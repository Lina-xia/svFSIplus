//xia

#ifndef COUPLE1D3D
#define COUPLE1D3D

class couple3D1D{
public:

    couple3D1D(int num){
        num_of_outlet = num;
        flow = new double[num_of_outlet];
        pressure = new double[num_of_outlet];
        // // test
        // for (int i = 0; i < num_of_outlet; i++){
        //     flow[i] = 2;
        // }
    }

    ~couple3D1D(){
        delete[] flow;
        delete[] pressure;
    }
 
    // 3D --> 1D
    static double* flow;
    // 1D --> 3D
    static double* pressure;

    static int num_of_outlet;
};
#endif