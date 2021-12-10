#pragma region SETUP
#include <iostream>
#include <fstream>
#include <queue>
#include <vector>
#include <string>
#include <deque>
#include <vector>
#include <assert.h>
using namespace std;
class Node;
typedef Node* Nodeptr;
#pragma endregion

#pragma region Modeling
typedef vector<float> fvector;

float dt = 0.001;
float current_time = 0;


void timeFlow()
{
    current_time += dt;
}

float getSignal() {
    if (current_time < 0)
        return 0;
    else return 1;
}

float rand_FloatRange(float a, float b)
{
    return ((b - a) * ((float)rand() / RAND_MAX)) + a;
}

class Node {
public:
    virtual float feedSignal(float input) = 0;
};

#pragma region GenericNodes

class GainNode : public Node {
public:
    float k;
    GainNode(float k) {
        this->k = k;
    }

    virtual float feedSignal(float input) {
        return input * k;
    }
};

class DiffNode : public Node {
public:
    virtual float feedSignal(float input) {
        float ret = (input - prev_input) / dt;
        prev_input = input;

        return ret;
    }

private:
    float prev_input = 0;
};

class DelayNode : public Node {
public:
    std::queue<float> _queue;

    float T = 0;
    DelayNode(float _T) : T(_T) {
    }

    virtual float feedSignal(float input) {

        _queue.push(input);

        float ret = 0;
        if (current_time > T)
        {
            ret = _queue.front();
            _queue.pop();
        }

        return ret;
    }
};

class IntNode : public Node {
public:

    IntNode(float k) {
        this->k = k;
    }

    virtual float feedSignal(float input) {

        int_value += 0.5 * dt * (input + prev_input) * k;
        prev_input = input;

        return int_value;
    }

private:
    float prev_input = 0;
    float int_value = 0;
    float k;
};

class SumNode : public Node {
public:

    SumNode(vector<float*> inputs, vector<int> signs)
        :signs(signs), second_inputs(inputs) {

        if (inputs.size() != signs.size())
            throw "Inputs and signs dimensions do not match in a summ node";

        for (int i = 0; i < signs.size(); i++)
            if (!(signs.at(i) == 1 || signs.at(i) == -1))
                throw "Invalid sign value in summ node";
    }

    virtual float feedSignal(float input) {
        float add = 0;

        for (int i = 0; i < second_inputs.size(); i++) {
            add += signs[i] * (*second_inputs[i]);
        }

        return input + add;
    }

private:
    vector<float*> second_inputs;
    vector<int> signs;
};

class TfNode : public Node
{
public:
    TfNode(fvector coefs)
        :coefs(coefs), prevVal(coefs.size()), first(1/coefs.at(0))
    {

    }

    virtual float feedSignal(float input) {

        float dirAndInputSum = input;

        for (int i = 1; i < prevVal.size(); i++)
            dirAndInputSum -= coefs.at(i) * prevVal.at(i);

        float lastDiriv = first * dirAndInputSum;

        prevVal[0] = lastDiriv;

        for (int i = 1; i < prevVal.size(); i++)
            prevVal[i] += prevVal.at(i - 1)*dt;

        return prevVal.at(prevVal.size() - 1);
    }

private:

    fvector coefs;
    fvector prevVal;
    float first;

};

class ForceNode : public Node
{
public:
    ForceNode(fvector coefs)
        :coefs(coefs), prevVal(coefs.size()), current(coefs.size())
    {

    }

    virtual float feedSignal(float input) {

        float out = coefs.at(coefs.size() - 1) * input;
        current[coefs.size() - 1] = input;

        for (int i = coefs.size() - 2; i >= 0; i--) {

            float diriv = (current.at(i+1) - prevVal.at(i+1)) / dt;
            out += coefs.at(i) * diriv;
            current[i] = diriv;

        }
        
        prevVal = current;
        return out;
    }

private:

    fvector coefs;
    fvector prevVal;
    fvector current;

};

class NoiseNode : public Node
{
public:
    NoiseNode(float m, float s)
        :mean(m), span(s/2)
    {

    }

    virtual float feedSignal(float input) {  
        return input + rand_FloatRange(mean - span, mean + span);
    }

private:

    float span;
    float mean;

};

class SaturationNode : public Node {
public:
    SaturationNode(float ins, float sut)
        :sut(sut),insens(ins)
    {

    }

    virtual float feedSignal(float input) {

        if (abs(input) < insens)
            return 0;
        else if (abs(input) >= sut)
            return sut;

        return input;
    }

private:

    float abs(float v) {
        return v < 0 ? -v : v;
    }
    float insens;
    float sut;

};

#pragma endregion


class NodeChain : public Node {
public:
    NodeChain(vector<Nodeptr> nodes) :nodes(nodes),len(nodes.size()) {

    }

    virtual float feedSignal(float input) {
        float ret = input;
        for (int i = 0; i < len; i++) {
            ret = nodes[i]->feedSignal(ret);
        }
        return ret;
    }

private:
    int len;
    vector<Nodeptr> nodes;
};

#pragma endregion

#pragma region Building



const string CIRCUT = "circut";
const string TF = "tf";
const string FORCE = "force";
const string GAIN = "gain";
const string INT = "int";
const string END = "endcircut";
const string NOISE = "noise";
const string SAT = "sat";

class Circut {
public:

    Circut(vector<Nodeptr> chain)
        : mainChain(chain)
    {

    }

    void Run(float maxT, ofstream& out ) {

        SumNode summ(vector<float*>{&outPoint}, vector<int>{-1});
        NodeChain chain(mainChain);

        while (current_time < maxT) {
            outPoint = chain.feedSignal(summ.feedSignal(getSignal()));
            out << current_time << '\t' << outPoint << endl;            
            timeFlow();
        }
    }

private:
    vector<Nodeptr> mainChain;
    float outPoint;
};

class Builder {

public:

    Builder(string filename)
        :input(filename) {
    }

    ~Builder() {
        input.close();
    }

    Circut Build() {
        string a = getNextWord();

        if(a == CIRCUT)
            BuildCircut();

        return Circut(mainChain);
    }

    void BuildCircut() {

        auto w = getNextWord();

        while (w != END) {
            if (w == TF)
                mainChain.push_back(BuildTf());
            else if (w == FORCE)
                mainChain.push_back(BuildForce());
            else if (w == GAIN)
                mainChain.push_back(BuildGain());
            else if (w == INT)
                mainChain.push_back(BuildInt());
            else if (w == SAT)
                mainChain.push_back(BuildSat());
            else if (w == NOISE)
                mainChain.push_back(BuildNoise());
            else {
                throw "";
            }

            w = getNextWord();
        }


    }

    Nodeptr BuildTf() {

        auto num = getNextWord();
        fvector coefs;

        while (num != ";") {
            coefs.push_back(atof(num.c_str()));
            num = getNextWord();
        }

        Nodeptr ret = new TfNode(coefs);
        return ret;
    }

    Nodeptr BuildForce() {

        auto num = getNextWord();
        fvector coefs;

        while (num != ";") {
            coefs.push_back(atof(num.c_str()));
            num = getNextWord();
        }

        Nodeptr ret = new ForceNode(coefs);
        return ret;
    }

    Nodeptr BuildGain() {

        auto num = getNextWord();
        fvector coefs;
        Nodeptr ret = new GainNode(atof(num.c_str()));
        return ret;
    }

    Nodeptr BuildInt() {

        auto num = getNextWord();
        fvector coefs;
        Nodeptr ret = new IntNode(atof(num.c_str()));
        return ret;
    }

    Nodeptr BuildSat() {

        auto ins = getNextWord();
        auto sat = getNextWord();

        Nodeptr ret = new SaturationNode(atof(ins.c_str()), atof(sat.c_str()));

        return ret;
    }

    Nodeptr BuildNoise() {

        auto mean = getNextWord();
        auto range = getNextWord();

        Nodeptr ret = new SaturationNode(atof(mean.c_str()), atof(range.c_str()));

        return ret;
    }

private:

    string getNextWord() { string a;  input >> a; return a; }

    ifstream input;
    vector<Nodeptr> mainChain;
    Nodeptr feedback;
};

#pragma endregion 

int main(int argc, char* argv[])
{
#pragma region INIT
    std::string out_filename;
    std::string in_filename;

    assert(argc == 3);

    in_filename = argv[1];
    out_filename = argv[2];

    int seed = 0;

#ifdef _DEBUG
    std::cout << "input  filename: " << in_filename << std::endl;
    std::cout << "output filename: " << out_filename << std::endl;
#endif // DEBUG

    
    for (int i = 0; i < in_filename.size(); i++)
        seed += static_cast<int>(in_filename.at(i));

    srand(seed);

    std::ofstream out(out_filename, 'w');


#pragma endregion

    Builder builder(in_filename);
    Circut c = builder.Build();
    c.Run(20, out);
  
    out.close();
}