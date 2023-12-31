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
rz(2.4581576) q[0];
sx q[0];
rz(12.0876) q[0];
rz(0.03102826) q[1];
sx q[1];
rz(-1.1614769) q[1];
sx q[1];
rz(-2.5002313) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.35533479) q[0];
sx q[0];
rz(-1.9541652) q[0];
sx q[0];
rz(-1.3134365) q[0];
rz(-pi) q[1];
rz(-1.9986721) q[2];
sx q[2];
rz(-1.2245721) q[2];
sx q[2];
rz(-1.8510173) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.76170834) q[1];
sx q[1];
rz(-2.0239081) q[1];
sx q[1];
rz(-0.74700991) q[1];
rz(-pi) q[2];
rz(0.62946837) q[3];
sx q[3];
rz(-1.9068204) q[3];
sx q[3];
rz(2.0365086) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(2.550094) q[2];
sx q[2];
rz(-1.9248362) q[2];
sx q[2];
rz(-0.47810289) q[2];
rz(1.6889307) q[3];
sx q[3];
rz(-2.0958459) q[3];
sx q[3];
rz(12/(7*pi)) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(pi/2) q[2];
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
rz(-0.61518413) q[1];
sx q[1];
rz(0.63308024) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.86242005) q[0];
sx q[0];
rz(-0.76403585) q[0];
sx q[0];
rz(0.11636244) q[0];
x q[1];
rz(2.5194397) q[2];
sx q[2];
rz(-2.2894147) q[2];
sx q[2];
rz(-2.6785786) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-2.9620348) q[1];
sx q[1];
rz(-1.9704559) q[1];
sx q[1];
rz(-2.7707997) q[1];
rz(-pi) q[2];
rz(-2.8301864) q[3];
sx q[3];
rz(-2.2283471) q[3];
sx q[3];
rz(-0.94535512) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-2.7818266) q[2];
sx q[2];
rz(-2.7414913) q[2];
sx q[2];
rz(0.11745545) q[2];
rz(-2.84058) q[3];
sx q[3];
rz(-1.8071226) q[3];
sx q[3];
rz(-1.7787748) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.85787073) q[0];
sx q[0];
rz(-2.4283333) q[0];
sx q[0];
rz(-0.088407956) q[0];
rz(1.1075426) q[1];
sx q[1];
rz(-2.2231367) q[1];
sx q[1];
rz(2.450768) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.847825) q[0];
sx q[0];
rz(-0.17058897) q[0];
sx q[0];
rz(-0.45844309) q[0];
rz(-pi) q[1];
rz(1.8059398) q[2];
sx q[2];
rz(-2.3217839) q[2];
sx q[2];
rz(-0.90607925) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(2.5882727) q[1];
sx q[1];
rz(-1.5447445) q[1];
sx q[1];
rz(-1.3887029) q[1];
rz(-pi) q[2];
x q[2];
rz(1.9224123) q[3];
sx q[3];
rz(-1.2120486) q[3];
sx q[3];
rz(-2.2946332) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-0.92418015) q[2];
sx q[2];
rz(-1.8841382) q[2];
sx q[2];
rz(-3.0351191) q[2];
rz(-1.4364093) q[3];
sx q[3];
rz(-0.36542106) q[3];
sx q[3];
rz(1.6586554) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0687662) q[0];
sx q[0];
rz(-2.9077353) q[0];
sx q[0];
rz(-2.6191214) q[0];
rz(2.8126295) q[1];
sx q[1];
rz(-1.4986228) q[1];
sx q[1];
rz(0.55363384) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5579917) q[0];
sx q[0];
rz(-0.86692536) q[0];
sx q[0];
rz(0.61917275) q[0];
x q[1];
rz(-0.76866863) q[2];
sx q[2];
rz(-0.93066356) q[2];
sx q[2];
rz(-1.4537571) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-1.1045038) q[1];
sx q[1];
rz(-1.4375086) q[1];
sx q[1];
rz(2.3549805) q[1];
rz(-pi) q[2];
x q[2];
rz(-3.1133075) q[3];
sx q[3];
rz(-2.3445498) q[3];
sx q[3];
rz(-2.459211) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.557495) q[2];
sx q[2];
rz(-0.9157052) q[2];
sx q[2];
rz(-2.6468357) q[2];
rz(-2.2385712) q[3];
sx q[3];
rz(-1.5938063) q[3];
sx q[3];
rz(1.2899227) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6506127) q[0];
sx q[0];
rz(-0.74478331) q[0];
sx q[0];
rz(-2.0565128) q[0];
rz(0.96013534) q[1];
sx q[1];
rz(-1.9624058) q[1];
sx q[1];
rz(0.18403149) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.40084307) q[0];
sx q[0];
rz(-1.6753734) q[0];
sx q[0];
rz(0.29086374) q[0];
rz(-pi) q[1];
x q[1];
rz(1.8475624) q[2];
sx q[2];
rz(-1.9230611) q[2];
sx q[2];
rz(0.89154746) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(2.9102238) q[1];
sx q[1];
rz(-2.0407045) q[1];
sx q[1];
rz(0.652657) q[1];
rz(-2.4155951) q[3];
sx q[3];
rz(-1.7222002) q[3];
sx q[3];
rz(-2.9263673) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(0.90157834) q[2];
sx q[2];
rz(-2.7928536) q[2];
sx q[2];
rz(-2.0007755) q[2];
rz(-2.7187637) q[3];
sx q[3];
rz(-1.4774277) q[3];
sx q[3];
rz(2.0146577) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.35247701) q[0];
sx q[0];
rz(-1.1081835) q[0];
sx q[0];
rz(-0.88678962) q[0];
rz(-2.9011762) q[1];
sx q[1];
rz(-1.0188894) q[1];
sx q[1];
rz(0.14850798) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7667023) q[0];
sx q[0];
rz(-2.3935211) q[0];
sx q[0];
rz(2.9038249) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.4609023) q[2];
sx q[2];
rz(-1.7981148) q[2];
sx q[2];
rz(2.0109039) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(3.0420694) q[1];
sx q[1];
rz(-1.6798786) q[1];
sx q[1];
rz(1.3101577) q[1];
rz(-pi) q[2];
rz(2.6753747) q[3];
sx q[3];
rz(-1.2122452) q[3];
sx q[3];
rz(2.8311604) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-0.85577661) q[2];
sx q[2];
rz(-0.4824051) q[2];
sx q[2];
rz(0.4883858) q[2];
rz(0.48940247) q[3];
sx q[3];
rz(-0.92178744) q[3];
sx q[3];
rz(1.2667806) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.69328904) q[0];
sx q[0];
rz(-1.332809) q[0];
sx q[0];
rz(0.81800246) q[0];
rz(1.8891634) q[1];
sx q[1];
rz(-2.7493582) q[1];
sx q[1];
rz(-2.0163527) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8669406) q[0];
sx q[0];
rz(-2.7526703) q[0];
sx q[0];
rz(1.2961943) q[0];
rz(2.0411885) q[2];
sx q[2];
rz(-0.33249582) q[2];
sx q[2];
rz(-1.5770797) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-0.69562558) q[1];
sx q[1];
rz(-1.274316) q[1];
sx q[1];
rz(-2.8545024) q[1];
rz(-pi) q[2];
rz(1.0057698) q[3];
sx q[3];
rz(-0.98494782) q[3];
sx q[3];
rz(0.98423959) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(2.0685048) q[2];
sx q[2];
rz(-2.1606074) q[2];
sx q[2];
rz(-2.4970064) q[2];
rz(-1.4792431) q[3];
sx q[3];
rz(-0.95696604) q[3];
sx q[3];
rz(-1.7361599) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1709764) q[0];
sx q[0];
rz(-3.0727486) q[0];
sx q[0];
rz(-1.6059426) q[0];
rz(-1.9203141) q[1];
sx q[1];
rz(-1.5131283) q[1];
sx q[1];
rz(-1.01064) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0897652) q[0];
sx q[0];
rz(-0.83501378) q[0];
sx q[0];
rz(-2.9249973) q[0];
rz(-pi) q[1];
x q[1];
rz(2.0660731) q[2];
sx q[2];
rz(-1.3066548) q[2];
sx q[2];
rz(1.2656982) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-1.1809363) q[1];
sx q[1];
rz(-1.8005383) q[1];
sx q[1];
rz(-2.9147663) q[1];
rz(0.40057064) q[3];
sx q[3];
rz(-1.6555602) q[3];
sx q[3];
rz(0.29464196) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-1.6179787) q[2];
sx q[2];
rz(-2.8221059) q[2];
sx q[2];
rz(-2.4475205) q[2];
rz(0.56898919) q[3];
sx q[3];
rz(-1.6945972) q[3];
sx q[3];
rz(-0.93769658) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0555608) q[0];
sx q[0];
rz(-1.9984351) q[0];
sx q[0];
rz(0.49474299) q[0];
rz(-2.5231979) q[1];
sx q[1];
rz(-1.6463552) q[1];
sx q[1];
rz(3.0659952) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9915174) q[0];
sx q[0];
rz(-0.65438327) q[0];
sx q[0];
rz(-2.0983216) q[0];
rz(-pi) q[1];
x q[1];
rz(3.0816239) q[2];
sx q[2];
rz(-1.7459933) q[2];
sx q[2];
rz(-0.8808459) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-0.52121431) q[1];
sx q[1];
rz(-1.673939) q[1];
sx q[1];
rz(2.4528273) q[1];
rz(-pi) q[2];
rz(-1.3887614) q[3];
sx q[3];
rz(-0.958003) q[3];
sx q[3];
rz(0.36422563) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(2.8081234) q[2];
sx q[2];
rz(-1.7227017) q[2];
sx q[2];
rz(1.3405651) q[2];
rz(-2.8373485) q[3];
sx q[3];
rz(-0.86383581) q[3];
sx q[3];
rz(0.58469599) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8463523) q[0];
sx q[0];
rz(-1.8363991) q[0];
sx q[0];
rz(-0.17679581) q[0];
rz(1.2416174) q[1];
sx q[1];
rz(-0.63443628) q[1];
sx q[1];
rz(-0.44100824) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2947185) q[0];
sx q[0];
rz(-0.28781578) q[0];
sx q[0];
rz(2.3527282) q[0];
rz(-pi) q[1];
rz(2.987864) q[2];
sx q[2];
rz(-1.4970461) q[2];
sx q[2];
rz(-1.7388625) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-2.0032138) q[1];
sx q[1];
rz(-0.56204501) q[1];
sx q[1];
rz(-2.2467062) q[1];
rz(-pi) q[2];
rz(1.8701843) q[3];
sx q[3];
rz(-0.96794879) q[3];
sx q[3];
rz(2.833948) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.56197721) q[2];
sx q[2];
rz(-2.5728971) q[2];
sx q[2];
rz(-2.8397172) q[2];
rz(0.89312303) q[3];
sx q[3];
rz(-1.8538657) q[3];
sx q[3];
rz(-2.9368403) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.72398913) q[0];
sx q[0];
rz(-1.2932734) q[0];
sx q[0];
rz(-1.4751157) q[0];
rz(0.026731116) q[1];
sx q[1];
rz(-1.4550799) q[1];
sx q[1];
rz(1.4310238) q[1];
rz(-1.0529636) q[2];
sx q[2];
rz(-0.79635194) q[2];
sx q[2];
rz(1.2780381) q[2];
rz(-2.5857153) q[3];
sx q[3];
rz(-2.0049958) q[3];
sx q[3];
rz(-0.47331664) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
