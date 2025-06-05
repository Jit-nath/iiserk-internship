#ifndef RECONSTRUCTOR_H
#define RECONSTRUCTOR_H
#include "include.h"
/*
                         ┌───────────────────────────┐
                         │        Alogorithm         │
                         └───────────────────────────┘
   ┌───────────────────────────────────────────────────────────────────────┐
   │ 1. Take the Atom name, position[x,y,z] and direction vector [x,y,z]   │
   │ 2. Take the reference graph which we want to re construct             │
   │ 3. read the reference graph                                           │
   │ 4. initialize a map {Atom id : {bond length, bond angle}}             │
   │ 5.
   │ 6.
   └───────────────────────────────────────────────────────────────────────┘
*/
class reconstructor
{
private:
public:
    static void PrintAlogrithm()
    {
        cout << "┌───────────────────────────┐" << endl;
        cout << "│        Alogorithm         │" << endl;
        cout << "└───────────────────────────┘" << endl;
        cout << endl;
        cout << "┌───────────────────────────────────────────────────────────────────────┐" << endl;
        cout << "│1. Take the Atom name, position[x,y,z] and direction vector [x,y,z]    │" << endl;
        cout << "│2. Take the reference graph which we want to re construct              │" << endl;
        cout << "│3. read the reference graph                                            │" << endl;
        cout << "│4. initialize a map {Atom id : {bond length, bond angle}}              │" << endl;
        cout << "│3.                                                                     │" << endl;
        cout << "└───────────────────────────────────────────────────────────────────────┘" << endl;
    }
};

#endif