OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[1],q[2],q[3];
rz(0.97857082) q[0];
sx q[0];
rz(-0.11712722) q[0];
sx q[0];
rz(0.22476619) q[0];
rz(2.6449142) q[1];
sx q[1];
rz(4.7159046) q[1];
sx q[1];
rz(11.386303) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.71043541) q[0];
sx q[0];
rz(-1.707516) q[0];
sx q[0];
rz(-1.4275292) q[0];
rz(0.013597537) q[2];
sx q[2];
rz(-0.44348479) q[2];
sx q[2];
rz(2.9421303) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(0.92946101) q[1];
sx q[1];
rz(-0.60570215) q[1];
sx q[1];
rz(-0.49537658) q[1];
rz(-1.7430844) q[3];
sx q[3];
rz(-2.7040561) q[3];
sx q[3];
rz(-0.019339081) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(2.4687389) q[2];
sx q[2];
rz(-1.3287975) q[2];
sx q[2];
rz(-1.9358181) q[2];
rz(-2.2852211) q[3];
sx q[3];
rz(-2.2090293) q[3];
sx q[3];
rz(0.87730733) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.67742753) q[0];
sx q[0];
rz(-2.7441315) q[0];
sx q[0];
rz(2.933266) q[0];
rz(-1.9147929) q[1];
sx q[1];
rz(-2.0715163) q[1];
sx q[1];
rz(2.5770381) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.78075942) q[0];
sx q[0];
rz(-1.6134904) q[0];
sx q[0];
rz(-1.4111931) q[0];
rz(2.4103569) q[2];
sx q[2];
rz(-1.1888767) q[2];
sx q[2];
rz(-2.3637091) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-2.1543252) q[1];
sx q[1];
rz(-0.72817737) q[1];
sx q[1];
rz(-0.2699234) q[1];
x q[2];
rz(-0.019314841) q[3];
sx q[3];
rz(-2.1105511) q[3];
sx q[3];
rz(3.0125953) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-0.63689256) q[2];
sx q[2];
rz(-1.8919287) q[2];
sx q[2];
rz(1.6553817) q[2];
rz(0.49556035) q[3];
sx q[3];
rz(-1.7335745) q[3];
sx q[3];
rz(1.0068309) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9946852) q[0];
sx q[0];
rz(-2.4201396) q[0];
sx q[0];
rz(1.7578693) q[0];
rz(1.9231298) q[1];
sx q[1];
rz(-1.0885025) q[1];
sx q[1];
rz(-2.7168435) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.95441636) q[0];
sx q[0];
rz(-1.5392051) q[0];
sx q[0];
rz(0.0626465) q[0];
rz(1.8239543) q[2];
sx q[2];
rz(-1.8033529) q[2];
sx q[2];
rz(-2.3311904) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(1.220881) q[1];
sx q[1];
rz(-1.6466116) q[1];
sx q[1];
rz(-0.37714173) q[1];
x q[2];
rz(-0.056428595) q[3];
sx q[3];
rz(-0.61739576) q[3];
sx q[3];
rz(-1.6625202) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(2.4212627) q[2];
sx q[2];
rz(-1.0676554) q[2];
sx q[2];
rz(0.024451582) q[2];
rz(-1.1484185) q[3];
sx q[3];
rz(-2.0423856) q[3];
sx q[3];
rz(-1.0816157) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9352683) q[0];
sx q[0];
rz(-0.22980389) q[0];
sx q[0];
rz(-0.21388737) q[0];
rz(0.81854406) q[1];
sx q[1];
rz(-1.7844424) q[1];
sx q[1];
rz(0.68665409) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.370433) q[0];
sx q[0];
rz(-2.4338182) q[0];
sx q[0];
rz(-2.4999655) q[0];
rz(-pi) q[1];
x q[1];
rz(1.0222203) q[2];
sx q[2];
rz(-1.8214263) q[2];
sx q[2];
rz(0.20285367) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-0.80685341) q[1];
sx q[1];
rz(-2.377393) q[1];
sx q[1];
rz(-3.0055226) q[1];
rz(-pi) q[2];
rz(-0.43709932) q[3];
sx q[3];
rz(-1.0954787) q[3];
sx q[3];
rz(0.10224414) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-2.116918) q[2];
sx q[2];
rz(-1.4093829) q[2];
sx q[2];
rz(0.68946687) q[2];
rz(1.8170478) q[3];
sx q[3];
rz(-1.9246512) q[3];
sx q[3];
rz(-0.57233468) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4631184) q[0];
sx q[0];
rz(-3.0480223) q[0];
sx q[0];
rz(2.1043188) q[0];
rz(2.4689238) q[1];
sx q[1];
rz(-1.9692407) q[1];
sx q[1];
rz(2.345828) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.138863) q[0];
sx q[0];
rz(-0.61289591) q[0];
sx q[0];
rz(-1.0010368) q[0];
rz(-pi) q[1];
rz(2.1433709) q[2];
sx q[2];
rz(-1.8279982) q[2];
sx q[2];
rz(-0.57359475) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(2.5368294) q[1];
sx q[1];
rz(-2.3336971) q[1];
sx q[1];
rz(1.1033586) q[1];
x q[2];
rz(-0.37973399) q[3];
sx q[3];
rz(-0.70443166) q[3];
sx q[3];
rz(-1.3317127) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-0.98046389) q[2];
sx q[2];
rz(-1.0995862) q[2];
sx q[2];
rz(0.83578342) q[2];
rz(-2.8849844) q[3];
sx q[3];
rz(-2.0308688) q[3];
sx q[3];
rz(-2.6512644) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4534054) q[0];
sx q[0];
rz(-1.9904933) q[0];
sx q[0];
rz(1.9292462) q[0];
rz(-2.1064099) q[1];
sx q[1];
rz(-1.6500708) q[1];
sx q[1];
rz(2.129668) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.69629151) q[0];
sx q[0];
rz(-0.20325771) q[0];
sx q[0];
rz(-2.1392034) q[0];
rz(-pi) q[1];
x q[1];
rz(1.5359319) q[2];
sx q[2];
rz(-2.8375585) q[2];
sx q[2];
rz(2.3140098) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-2.1862667) q[1];
sx q[1];
rz(-1.9791934) q[1];
sx q[1];
rz(-0.46640654) q[1];
rz(-pi) q[2];
rz(-1.8406781) q[3];
sx q[3];
rz(-1.9414475) q[3];
sx q[3];
rz(-2.5605367) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-2.6637491) q[2];
sx q[2];
rz(-2.9426844) q[2];
sx q[2];
rz(0.94856962) q[2];
rz(2.6712724) q[3];
sx q[3];
rz(-2.0430653) q[3];
sx q[3];
rz(-1.8868014) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.54742852) q[0];
sx q[0];
rz(-2.1305888) q[0];
sx q[0];
rz(1.2642566) q[0];
rz(-2.3612379) q[1];
sx q[1];
rz(-0.64096132) q[1];
sx q[1];
rz(2.9497214) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.22671285) q[0];
sx q[0];
rz(-0.60189819) q[0];
sx q[0];
rz(-2.4826217) q[0];
rz(-0.35991798) q[2];
sx q[2];
rz(-1.214817) q[2];
sx q[2];
rz(2.6448665) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(2.807232) q[1];
sx q[1];
rz(-1.8340115) q[1];
sx q[1];
rz(0.46137793) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.97080462) q[3];
sx q[3];
rz(-1.8626584) q[3];
sx q[3];
rz(-1.6455022) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-0.65036217) q[2];
sx q[2];
rz(-1.0543084) q[2];
sx q[2];
rz(0.80797705) q[2];
rz(-2.2982277) q[3];
sx q[3];
rz(-1.9873514) q[3];
sx q[3];
rz(-1.556506) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4845487) q[0];
sx q[0];
rz(-0.024024155) q[0];
sx q[0];
rz(2.3727544) q[0];
rz(2.8688042) q[1];
sx q[1];
rz(-1.6285248) q[1];
sx q[1];
rz(1.3974894) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.63015712) q[0];
sx q[0];
rz(-0.4841329) q[0];
sx q[0];
rz(0.70229437) q[0];
x q[1];
rz(-2.0968123) q[2];
sx q[2];
rz(-0.57990852) q[2];
sx q[2];
rz(-1.8749664) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(0.2915512) q[1];
sx q[1];
rz(-1.631389) q[1];
sx q[1];
rz(-1.4539976) q[1];
x q[2];
rz(0.66889735) q[3];
sx q[3];
rz(-1.9188768) q[3];
sx q[3];
rz(-2.477248) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-2.2046795) q[2];
sx q[2];
rz(-1.0656554) q[2];
sx q[2];
rz(1.6395456) q[2];
rz(1.8223193) q[3];
sx q[3];
rz(-2.2917512) q[3];
sx q[3];
rz(1.5892861) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.69916344) q[0];
sx q[0];
rz(-0.84446877) q[0];
sx q[0];
rz(-0.97802877) q[0];
rz(-2.6507822) q[1];
sx q[1];
rz(-1.4421137) q[1];
sx q[1];
rz(1.6455654) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.42616823) q[0];
sx q[0];
rz(-1.4472571) q[0];
sx q[0];
rz(-2.9036595) q[0];
rz(1.3991576) q[2];
sx q[2];
rz(-1.1632068) q[2];
sx q[2];
rz(-1.379058) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.2435808) q[1];
sx q[1];
rz(-1.812938) q[1];
sx q[1];
rz(-2.8088074) q[1];
x q[2];
rz(1.6469025) q[3];
sx q[3];
rz(-1.0400912) q[3];
sx q[3];
rz(0.87065856) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(0.28828037) q[2];
sx q[2];
rz(-2.2744982) q[2];
sx q[2];
rz(-1.3221928) q[2];
rz(0.36618048) q[3];
sx q[3];
rz(-1.7145232) q[3];
sx q[3];
rz(1.1314932) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.91728297) q[0];
sx q[0];
rz(-1.5320822) q[0];
sx q[0];
rz(-3.068058) q[0];
rz(-1.3167943) q[1];
sx q[1];
rz(-1.2578745) q[1];
sx q[1];
rz(-1.8427461) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6599346) q[0];
sx q[0];
rz(-1.3941688) q[0];
sx q[0];
rz(0.28089471) q[0];
rz(-pi) q[1];
rz(-1.3234947) q[2];
sx q[2];
rz(-2.6404233) q[2];
sx q[2];
rz(2.4275818) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-1.6886814) q[1];
sx q[1];
rz(-1.9454524) q[1];
sx q[1];
rz(-0.25325798) q[1];
x q[2];
rz(0.73963005) q[3];
sx q[3];
rz(-1.6508023) q[3];
sx q[3];
rz(-2.7440939) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(1.9650044) q[2];
sx q[2];
rz(-2.1829288) q[2];
sx q[2];
rz(2.5625572) q[2];
rz(-1.6252919) q[3];
sx q[3];
rz(-0.54280353) q[3];
sx q[3];
rz(-1.6453086) q[3];
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
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8662921) q[0];
sx q[0];
rz(-1.0466812) q[0];
sx q[0];
rz(-0.064591797) q[0];
rz(-2.2304089) q[1];
sx q[1];
rz(-1.8235527) q[1];
sx q[1];
rz(1.8240737) q[1];
rz(-0.42413128) q[2];
sx q[2];
rz(-1.567622) q[2];
sx q[2];
rz(1.4040053) q[2];
rz(1.961962) q[3];
sx q[3];
rz(-2.4868271) q[3];
sx q[3];
rz(0.42533608) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
