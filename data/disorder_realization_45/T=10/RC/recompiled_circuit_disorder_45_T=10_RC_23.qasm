OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-1.835445) q[0];
sx q[0];
rz(-0.68343502) q[0];
sx q[0];
rz(0.47877065) q[0];
rz(-3.1105644) q[1];
sx q[1];
rz(-1.9801158) q[1];
sx q[1];
rz(2.5002313) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7862579) q[0];
sx q[0];
rz(-1.9541652) q[0];
sx q[0];
rz(-1.3134365) q[0];
rz(-1.1429206) q[2];
sx q[2];
rz(-1.2245721) q[2];
sx q[2];
rz(-1.2905754) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(0.42395017) q[1];
sx q[1];
rz(-0.91361928) q[1];
sx q[1];
rz(0.98492019) q[1];
x q[2];
rz(0.53545973) q[3];
sx q[3];
rz(-0.70264953) q[3];
sx q[3];
rz(3.1010166) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-2.550094) q[2];
sx q[2];
rz(-1.9248362) q[2];
sx q[2];
rz(-2.6634898) q[2];
rz(-1.6889307) q[3];
sx q[3];
rz(-1.0457467) q[3];
sx q[3];
rz(-2.5959192) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.96145445) q[0];
sx q[0];
rz(-0.77471662) q[0];
sx q[0];
rz(2.1226728) q[0];
rz(-1.6628751) q[1];
sx q[1];
rz(-2.5264085) q[1];
sx q[1];
rz(2.5085124) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.7019254) q[0];
sx q[0];
rz(-0.81322008) q[0];
sx q[0];
rz(1.6815835) q[0];
x q[1];
rz(2.1585629) q[2];
sx q[2];
rz(-2.2289742) q[2];
sx q[2];
rz(-0.36511974) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-0.96453297) q[1];
sx q[1];
rz(-2.6033083) q[1];
sx q[1];
rz(2.2798377) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.8301864) q[3];
sx q[3];
rz(-2.2283471) q[3];
sx q[3];
rz(-0.94535512) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-0.35976609) q[2];
sx q[2];
rz(-0.40010139) q[2];
sx q[2];
rz(-3.0241372) q[2];
rz(-2.84058) q[3];
sx q[3];
rz(-1.3344701) q[3];
sx q[3];
rz(1.7787748) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.85787073) q[0];
sx q[0];
rz(-2.4283333) q[0];
sx q[0];
rz(0.088407956) q[0];
rz(-1.1075426) q[1];
sx q[1];
rz(-2.2231367) q[1];
sx q[1];
rz(0.69082469) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.847825) q[0];
sx q[0];
rz(-2.9710037) q[0];
sx q[0];
rz(-2.6831496) q[0];
rz(-pi) q[1];
rz(-2.3766741) q[2];
sx q[2];
rz(-1.3996482) q[2];
sx q[2];
rz(-2.314832) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-1.1580116) q[1];
sx q[1];
rz(-2.9576655) q[1];
sx q[1];
rz(-1.4278825) q[1];
rz(-pi) q[2];
x q[2];
rz(2.7615943) q[3];
sx q[3];
rz(-1.899154) q[3];
sx q[3];
rz(0.85193714) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.2174125) q[2];
sx q[2];
rz(-1.2574544) q[2];
sx q[2];
rz(0.10647354) q[2];
rz(-1.7051833) q[3];
sx q[3];
rz(-0.36542106) q[3];
sx q[3];
rz(1.4829372) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0728264) q[0];
sx q[0];
rz(-2.9077353) q[0];
sx q[0];
rz(-0.52247125) q[0];
rz(-0.32896313) q[1];
sx q[1];
rz(-1.4986228) q[1];
sx q[1];
rz(-2.5879588) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.41933665) q[0];
sx q[0];
rz(-1.1124848) q[0];
sx q[0];
rz(-0.76461794) q[0];
x q[1];
rz(-0.81981084) q[2];
sx q[2];
rz(-2.1852583) q[2];
sx q[2];
rz(-2.7053506) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.33448996) q[1];
sx q[1];
rz(-2.346171) q[1];
sx q[1];
rz(2.9544178) q[1];
x q[2];
rz(0.028285154) q[3];
sx q[3];
rz(-0.79704282) q[3];
sx q[3];
rz(-0.68238168) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-0.557495) q[2];
sx q[2];
rz(-0.9157052) q[2];
sx q[2];
rz(-2.6468357) q[2];
rz(2.2385712) q[3];
sx q[3];
rz(-1.5938063) q[3];
sx q[3];
rz(-1.2899227) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.49098) q[0];
sx q[0];
rz(-2.3968093) q[0];
sx q[0];
rz(2.0565128) q[0];
rz(2.1814573) q[1];
sx q[1];
rz(-1.1791869) q[1];
sx q[1];
rz(0.18403149) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.40084307) q[0];
sx q[0];
rz(-1.4662192) q[0];
sx q[0];
rz(-2.8507289) q[0];
rz(-pi) q[1];
rz(-0.63921914) q[2];
sx q[2];
rz(-2.6972065) q[2];
sx q[2];
rz(1.5803312) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(2.3369493) q[1];
sx q[1];
rz(-0.78362432) q[1];
sx q[1];
rz(0.6964535) q[1];
rz(2.4155951) q[3];
sx q[3];
rz(-1.4193924) q[3];
sx q[3];
rz(0.21522537) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-2.2400143) q[2];
sx q[2];
rz(-0.34873909) q[2];
sx q[2];
rz(-1.1408172) q[2];
rz(0.42282894) q[3];
sx q[3];
rz(-1.6641649) q[3];
sx q[3];
rz(1.1269349) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.35247701) q[0];
sx q[0];
rz(-1.1081835) q[0];
sx q[0];
rz(-0.88678962) q[0];
rz(-0.24041644) q[1];
sx q[1];
rz(-2.1227032) q[1];
sx q[1];
rz(-2.9930847) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1214971) q[0];
sx q[0];
rz(-1.4098865) q[0];
sx q[0];
rz(0.73385977) q[0];
rz(1.4609023) q[2];
sx q[2];
rz(-1.3434778) q[2];
sx q[2];
rz(2.0109039) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.099523274) q[1];
sx q[1];
rz(-1.4617141) q[1];
sx q[1];
rz(1.3101577) q[1];
rz(-pi) q[2];
rz(0.69494436) q[3];
sx q[3];
rz(-2.5616025) q[3];
sx q[3];
rz(0.65141962) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(2.285816) q[2];
sx q[2];
rz(-2.6591876) q[2];
sx q[2];
rz(-0.4883858) q[2];
rz(-2.6521902) q[3];
sx q[3];
rz(-0.92178744) q[3];
sx q[3];
rz(-1.874812) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.69328904) q[0];
sx q[0];
rz(-1.332809) q[0];
sx q[0];
rz(0.81800246) q[0];
rz(-1.8891634) q[1];
sx q[1];
rz(-2.7493582) q[1];
sx q[1];
rz(2.0163527) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.274652) q[0];
sx q[0];
rz(-2.7526703) q[0];
sx q[0];
rz(1.8453983) q[0];
rz(-0.15525012) q[2];
sx q[2];
rz(-1.2755738) q[2];
sx q[2];
rz(-2.070602) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.69562558) q[1];
sx q[1];
rz(-1.274316) q[1];
sx q[1];
rz(-0.28709025) q[1];
x q[2];
rz(0.66594395) q[3];
sx q[3];
rz(-1.1083372) q[3];
sx q[3];
rz(2.8921814) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.0730878) q[2];
sx q[2];
rz(-0.98098522) q[2];
sx q[2];
rz(-2.4970064) q[2];
rz(-1.4792431) q[3];
sx q[3];
rz(-2.1846266) q[3];
sx q[3];
rz(-1.4054327) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.97061625) q[0];
sx q[0];
rz(-0.068844065) q[0];
sx q[0];
rz(-1.53565) q[0];
rz(1.9203141) q[1];
sx q[1];
rz(-1.6284643) q[1];
sx q[1];
rz(2.1309526) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0897652) q[0];
sx q[0];
rz(-0.83501378) q[0];
sx q[0];
rz(2.9249973) q[0];
x q[1];
rz(2.0881565) q[2];
sx q[2];
rz(-0.55609716) q[2];
sx q[2];
rz(-2.9964787) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(1.9606564) q[1];
sx q[1];
rz(-1.3410543) q[1];
sx q[1];
rz(-0.22682637) q[1];
rz(-pi) q[2];
rz(0.2145433) q[3];
sx q[3];
rz(-0.40896591) q[3];
sx q[3];
rz(1.4734801) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-1.5236139) q[2];
sx q[2];
rz(-0.31948677) q[2];
sx q[2];
rz(-2.4475205) q[2];
rz(2.5726035) q[3];
sx q[3];
rz(-1.4469955) q[3];
sx q[3];
rz(2.2038961) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.086031832) q[0];
sx q[0];
rz(-1.9984351) q[0];
sx q[0];
rz(2.6468497) q[0];
rz(-2.5231979) q[1];
sx q[1];
rz(-1.6463552) q[1];
sx q[1];
rz(3.0659952) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.15007524) q[0];
sx q[0];
rz(-0.65438327) q[0];
sx q[0];
rz(-1.0432711) q[0];
x q[1];
rz(-1.2443301) q[2];
sx q[2];
rz(-0.18507659) q[2];
sx q[2];
rz(-1.9290123) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.52121431) q[1];
sx q[1];
rz(-1.673939) q[1];
sx q[1];
rz(0.68876536) q[1];
x q[2];
rz(-1.7528312) q[3];
sx q[3];
rz(-2.1835897) q[3];
sx q[3];
rz(0.36422563) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-0.33346924) q[2];
sx q[2];
rz(-1.418891) q[2];
sx q[2];
rz(1.8010275) q[2];
rz(0.30424413) q[3];
sx q[3];
rz(-0.86383581) q[3];
sx q[3];
rz(-2.5568967) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8463523) q[0];
sx q[0];
rz(-1.3051935) q[0];
sx q[0];
rz(2.9647968) q[0];
rz(1.8999752) q[1];
sx q[1];
rz(-0.63443628) q[1];
sx q[1];
rz(0.44100824) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.037017578) q[0];
sx q[0];
rz(-1.772176) q[0];
sx q[0];
rz(-1.3637533) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.6454234) q[2];
sx q[2];
rz(-1.4174889) q[2];
sx q[2];
rz(0.15664936) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-0.16377201) q[1];
sx q[1];
rz(-1.2308916) q[1];
sx q[1];
rz(2.0274859) q[1];
rz(-pi) q[2];
rz(-2.7367758) q[3];
sx q[3];
rz(-2.4768618) q[3];
sx q[3];
rz(2.3354195) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(2.5796154) q[2];
sx q[2];
rz(-2.5728971) q[2];
sx q[2];
rz(2.8397172) q[2];
rz(-2.2484696) q[3];
sx q[3];
rz(-1.2877269) q[3];
sx q[3];
rz(-0.20475234) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4176035) q[0];
sx q[0];
rz(-1.2932734) q[0];
sx q[0];
rz(-1.4751157) q[0];
rz(-3.1148615) q[1];
sx q[1];
rz(-1.4550799) q[1];
sx q[1];
rz(1.4310238) q[1];
rz(-1.0529636) q[2];
sx q[2];
rz(-0.79635194) q[2];
sx q[2];
rz(1.2780381) q[2];
rz(2.0704913) q[3];
sx q[3];
rz(-1.0715967) q[3];
sx q[3];
rz(-1.788492) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
