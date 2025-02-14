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
rz(0.12953144) q[0];
sx q[0];
rz(3.01053) q[0];
sx q[0];
rz(13.170903) q[0];
rz(-0.52855748) q[1];
sx q[1];
rz(-0.25783917) q[1];
sx q[1];
rz(-1.2300307) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.226868) q[0];
sx q[0];
rz(-0.42166963) q[0];
sx q[0];
rz(-0.64897169) q[0];
rz(2.8240439) q[2];
sx q[2];
rz(-2.5032836) q[2];
sx q[2];
rz(-2.2090863) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.22495077) q[1];
sx q[1];
rz(-1.0826546) q[1];
sx q[1];
rz(0.78190885) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.3377519) q[3];
sx q[3];
rz(-2.6342158) q[3];
sx q[3];
rz(-0.22358596) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-1.0543542) q[2];
sx q[2];
rz(-1.9638991) q[2];
sx q[2];
rz(-0.13620201) q[2];
rz(-0.44998351) q[3];
sx q[3];
rz(-2.4583702) q[3];
sx q[3];
rz(2.3581678) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0403274) q[0];
sx q[0];
rz(-0.7985006) q[0];
sx q[0];
rz(-0.26327565) q[0];
rz(-2.1030262) q[1];
sx q[1];
rz(-2.7949605) q[1];
sx q[1];
rz(2.3668049) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2939071) q[0];
sx q[0];
rz(-0.65252248) q[0];
sx q[0];
rz(0.34282617) q[0];
rz(-pi) q[1];
rz(-0.75157849) q[2];
sx q[2];
rz(-2.2937991) q[2];
sx q[2];
rz(2.4931049) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-2.0940184) q[1];
sx q[1];
rz(-1.9500004) q[1];
sx q[1];
rz(-3.1017041) q[1];
x q[2];
rz(2.5790096) q[3];
sx q[3];
rz(-2.354971) q[3];
sx q[3];
rz(-0.41401573) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.420555) q[2];
sx q[2];
rz(-1.5169531) q[2];
sx q[2];
rz(0.59845412) q[2];
rz(1.7789486) q[3];
sx q[3];
rz(-2.2612031) q[3];
sx q[3];
rz(-2.8708598) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8552928) q[0];
sx q[0];
rz(-2.6787651) q[0];
sx q[0];
rz(1.5077952) q[0];
rz(-0.15448054) q[1];
sx q[1];
rz(-1.721761) q[1];
sx q[1];
rz(-2.3522164) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.75784439) q[0];
sx q[0];
rz(-2.0998635) q[0];
sx q[0];
rz(1.4704513) q[0];
x q[1];
rz(-0.28114281) q[2];
sx q[2];
rz(-1.707555) q[2];
sx q[2];
rz(-0.006055486) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(2.912108) q[1];
sx q[1];
rz(-2.4463725) q[1];
sx q[1];
rz(2.0575614) q[1];
rz(-pi) q[2];
x q[2];
rz(2.0752729) q[3];
sx q[3];
rz(-2.3036727) q[3];
sx q[3];
rz(1.63597) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(0.35884759) q[2];
sx q[2];
rz(-2.2300356) q[2];
sx q[2];
rz(-1.6645128) q[2];
rz(1.9278256) q[3];
sx q[3];
rz(-1.8002847) q[3];
sx q[3];
rz(1.7267797) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3690935) q[0];
sx q[0];
rz(-1.0558244) q[0];
sx q[0];
rz(-0.61071998) q[0];
rz(0.67131132) q[1];
sx q[1];
rz(-1.5076312) q[1];
sx q[1];
rz(-1.3549365) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8313361) q[0];
sx q[0];
rz(-0.68643565) q[0];
sx q[0];
rz(-2.3879334) q[0];
x q[1];
rz(-1.942039) q[2];
sx q[2];
rz(-1.4535731) q[2];
sx q[2];
rz(3.0407112) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(0.23933168) q[1];
sx q[1];
rz(-0.25333764) q[1];
sx q[1];
rz(1.317666) q[1];
rz(-pi) q[2];
x q[2];
rz(0.98396222) q[3];
sx q[3];
rz(-1.6500743) q[3];
sx q[3];
rz(-0.51494277) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-1.9352202) q[2];
sx q[2];
rz(-0.78038961) q[2];
sx q[2];
rz(2.8727403) q[2];
rz(2.3937461) q[3];
sx q[3];
rz(-0.81955376) q[3];
sx q[3];
rz(-1.4486754) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.95626962) q[0];
sx q[0];
rz(-1.3893501) q[0];
sx q[0];
rz(2.4638033) q[0];
rz(2.9970844) q[1];
sx q[1];
rz(-0.78737193) q[1];
sx q[1];
rz(-1.6544624) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5337473) q[0];
sx q[0];
rz(-0.13838875) q[0];
sx q[0];
rz(1.3800623) q[0];
rz(0.057737902) q[2];
sx q[2];
rz(-1.9721834) q[2];
sx q[2];
rz(2.558208) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-2.2320391) q[1];
sx q[1];
rz(-2.0318446) q[1];
sx q[1];
rz(-0.98566815) q[1];
rz(1.1950567) q[3];
sx q[3];
rz(-1.5259966) q[3];
sx q[3];
rz(-0.57023772) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(0.80374485) q[2];
sx q[2];
rz(-0.2636815) q[2];
sx q[2];
rz(2.7860876) q[2];
rz(-1.0454987) q[3];
sx q[3];
rz(-1.9020566) q[3];
sx q[3];
rz(-2.579328) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1657555) q[0];
sx q[0];
rz(-2.5582357) q[0];
sx q[0];
rz(-3.052886) q[0];
rz(-1.6070131) q[1];
sx q[1];
rz(-1.8314223) q[1];
sx q[1];
rz(-3.065899) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.59654419) q[0];
sx q[0];
rz(-0.91191429) q[0];
sx q[0];
rz(1.8152899) q[0];
rz(-pi) q[1];
rz(1.6275109) q[2];
sx q[2];
rz(-0.73243388) q[2];
sx q[2];
rz(2.5947941) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-2.9348889) q[1];
sx q[1];
rz(-2.5674501) q[1];
sx q[1];
rz(-1.8341792) q[1];
rz(-2.663732) q[3];
sx q[3];
rz(-0.63242542) q[3];
sx q[3];
rz(-0.53148182) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(2.3169516) q[2];
sx q[2];
rz(-1.3836766) q[2];
sx q[2];
rz(-1.0392044) q[2];
rz(-0.21229395) q[3];
sx q[3];
rz(-1.1024691) q[3];
sx q[3];
rz(1.4958517) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.067326389) q[0];
sx q[0];
rz(-2.3889611) q[0];
sx q[0];
rz(2.0972032) q[0];
rz(2.3692865) q[1];
sx q[1];
rz(-2.4817395) q[1];
sx q[1];
rz(-2.4078802) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.70673101) q[0];
sx q[0];
rz(-1.0821663) q[0];
sx q[0];
rz(2.3948993) q[0];
x q[1];
rz(2.9915221) q[2];
sx q[2];
rz(-1.2037306) q[2];
sx q[2];
rz(2.31524) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-0.83625863) q[1];
sx q[1];
rz(-1.6871258) q[1];
sx q[1];
rz(-0.29246026) q[1];
rz(1.4513408) q[3];
sx q[3];
rz(-2.4032695) q[3];
sx q[3];
rz(-2.8286162) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(2.5250728) q[2];
sx q[2];
rz(-2.6313582) q[2];
sx q[2];
rz(-0.60866848) q[2];
rz(1.0673374) q[3];
sx q[3];
rz(-2.267024) q[3];
sx q[3];
rz(2.5909891) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6390425) q[0];
sx q[0];
rz(-2.783343) q[0];
sx q[0];
rz(-1.6453561) q[0];
rz(-1.3642338) q[1];
sx q[1];
rz(-0.61449209) q[1];
sx q[1];
rz(-0.38633698) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5674114) q[0];
sx q[0];
rz(-1.8339388) q[0];
sx q[0];
rz(2.4469923) q[0];
rz(-pi) q[1];
rz(-1.3851829) q[2];
sx q[2];
rz(-1.6395373) q[2];
sx q[2];
rz(0.5663213) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(2.1586502) q[1];
sx q[1];
rz(-1.6931931) q[1];
sx q[1];
rz(2.8084408) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.59054116) q[3];
sx q[3];
rz(-1.6034551) q[3];
sx q[3];
rz(1.5368568) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(2.6284457) q[2];
sx q[2];
rz(-1.7690423) q[2];
sx q[2];
rz(2.9885542) q[2];
rz(0.89093527) q[3];
sx q[3];
rz(-0.95518437) q[3];
sx q[3];
rz(2.4002767) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1560169) q[0];
sx q[0];
rz(-0.22933904) q[0];
sx q[0];
rz(-0.34737059) q[0];
rz(0.037847606) q[1];
sx q[1];
rz(-1.9719351) q[1];
sx q[1];
rz(-2.6191424) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.98451) q[0];
sx q[0];
rz(-1.1584131) q[0];
sx q[0];
rz(-0.097481485) q[0];
rz(-pi) q[1];
x q[1];
rz(0.36549296) q[2];
sx q[2];
rz(-1.5710982) q[2];
sx q[2];
rz(-3.0916758) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(0.47679893) q[1];
sx q[1];
rz(-2.0719686) q[1];
sx q[1];
rz(-1.0887926) q[1];
rz(-0.40710514) q[3];
sx q[3];
rz(-2.0927755) q[3];
sx q[3];
rz(0.84195053) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-0.6309506) q[2];
sx q[2];
rz(-0.46322552) q[2];
sx q[2];
rz(1.9412712) q[2];
rz(0.55780324) q[3];
sx q[3];
rz(-1.5188981) q[3];
sx q[3];
rz(-2.7571078) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2931622) q[0];
sx q[0];
rz(-0.98386216) q[0];
sx q[0];
rz(-1.7278607) q[0];
rz(-3.0005455) q[1];
sx q[1];
rz(-1.3755362) q[1];
sx q[1];
rz(0.65473762) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.075874627) q[0];
sx q[0];
rz(-1.362934) q[0];
sx q[0];
rz(-1.8130568) q[0];
rz(-pi) q[1];
rz(1.0626002) q[2];
sx q[2];
rz(-2.104665) q[2];
sx q[2];
rz(-1.2956308) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(2.4550528) q[1];
sx q[1];
rz(-1.3994263) q[1];
sx q[1];
rz(1.8167348) q[1];
x q[2];
rz(1.039417) q[3];
sx q[3];
rz(-2.1221042) q[3];
sx q[3];
rz(1.2599864) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-0.35141382) q[2];
sx q[2];
rz(-1.979579) q[2];
sx q[2];
rz(-0.72511017) q[2];
rz(0.890598) q[3];
sx q[3];
rz(-2.2913439) q[3];
sx q[3];
rz(-2.5543673) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.17726041) q[0];
sx q[0];
rz(-1.6163419) q[0];
sx q[0];
rz(1.1905715) q[0];
rz(1.1217077) q[1];
sx q[1];
rz(-1.5737166) q[1];
sx q[1];
rz(2.166688) q[1];
rz(-1.2091985) q[2];
sx q[2];
rz(-0.93930106) q[2];
sx q[2];
rz(1.3472547) q[2];
rz(2.4287281) q[3];
sx q[3];
rz(-1.6920857) q[3];
sx q[3];
rz(-1.6209775) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
