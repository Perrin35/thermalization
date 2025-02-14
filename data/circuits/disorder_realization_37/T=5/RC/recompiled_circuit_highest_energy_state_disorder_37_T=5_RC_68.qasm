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
rz(-0.77684075) q[0];
sx q[0];
rz(2.2651894) q[0];
sx q[0];
rz(9.6776008) q[0];
rz(1.1480992) q[1];
sx q[1];
rz(4.0568772) q[1];
sx q[1];
rz(8.6203909) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0067449) q[0];
sx q[0];
rz(-2.6383556) q[0];
sx q[0];
rz(-0.87849599) q[0];
rz(-pi) q[1];
x q[1];
rz(0.47801272) q[2];
sx q[2];
rz(-0.4768663) q[2];
sx q[2];
rz(0.074772686) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(1.8411257) q[1];
sx q[1];
rz(-1.2004832) q[1];
sx q[1];
rz(-1.9410067) q[1];
rz(-pi) q[2];
x q[2];
rz(0.34396307) q[3];
sx q[3];
rz(-1.4620145) q[3];
sx q[3];
rz(3.1358842) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-3.048451) q[2];
sx q[2];
rz(-0.96258771) q[2];
sx q[2];
rz(0.87240458) q[2];
rz(-2.8060272) q[3];
sx q[3];
rz(-1.395697) q[3];
sx q[3];
rz(-0.083757639) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
rz(-pi) q[3];
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
rz(0.63808477) q[0];
sx q[0];
rz(-2.8758949) q[0];
sx q[0];
rz(2.9124394) q[0];
rz(-2.9511662) q[1];
sx q[1];
rz(-1.3399905) q[1];
sx q[1];
rz(0.71358877) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.021727) q[0];
sx q[0];
rz(-1.83868) q[0];
sx q[0];
rz(-2.2557206) q[0];
x q[1];
rz(-0.53411463) q[2];
sx q[2];
rz(-2.2306109) q[2];
sx q[2];
rz(2.5995863) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-1.16878) q[1];
sx q[1];
rz(-2.3832364) q[1];
sx q[1];
rz(0.43557628) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.3768164) q[3];
sx q[3];
rz(-1.4524609) q[3];
sx q[3];
rz(-1.3142933) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(2.4552292) q[2];
sx q[2];
rz(-2.8296622) q[2];
sx q[2];
rz(2.6658106) q[2];
rz(-1.0698211) q[3];
sx q[3];
rz(-2.1333623) q[3];
sx q[3];
rz(-0.76655918) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0723202) q[0];
sx q[0];
rz(-1.197149) q[0];
sx q[0];
rz(1.9092165) q[0];
rz(0.45066372) q[1];
sx q[1];
rz(-1.9602937) q[1];
sx q[1];
rz(2.252069) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.33143932) q[0];
sx q[0];
rz(-0.74140775) q[0];
sx q[0];
rz(-0.62192328) q[0];
rz(-pi) q[1];
x q[1];
rz(0.53307791) q[2];
sx q[2];
rz(-2.0642082) q[2];
sx q[2];
rz(2.5149038) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-1.3865859) q[1];
sx q[1];
rz(-0.68803794) q[1];
sx q[1];
rz(-1.2351456) q[1];
rz(2.398185) q[3];
sx q[3];
rz(-2.6498859) q[3];
sx q[3];
rz(-1.7984901) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(1.4543317) q[2];
sx q[2];
rz(-0.4250409) q[2];
sx q[2];
rz(-1.8827776) q[2];
rz(1.2339633) q[3];
sx q[3];
rz(-1.1326658) q[3];
sx q[3];
rz(-0.0080571938) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.71037978) q[0];
sx q[0];
rz(-2.2358535) q[0];
sx q[0];
rz(1.4696962) q[0];
rz(-1.1471033) q[1];
sx q[1];
rz(-2.1397782) q[1];
sx q[1];
rz(-0.9009487) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9757864) q[0];
sx q[0];
rz(-0.32484178) q[0];
sx q[0];
rz(-2.8715517) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.5261648) q[2];
sx q[2];
rz(-1.0503029) q[2];
sx q[2];
rz(2.8857638) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(1.683953) q[1];
sx q[1];
rz(-1.6648653) q[1];
sx q[1];
rz(2.4165618) q[1];
rz(-2.4086558) q[3];
sx q[3];
rz(-1.4519435) q[3];
sx q[3];
rz(-0.87866966) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-2.6441696) q[2];
sx q[2];
rz(-1.7449417) q[2];
sx q[2];
rz(0.72745848) q[2];
rz(-2.4265031) q[3];
sx q[3];
rz(-2.6810724) q[3];
sx q[3];
rz(-2.6434744) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.46399507) q[0];
sx q[0];
rz(-1.3901187) q[0];
sx q[0];
rz(0.736262) q[0];
rz(2.5212506) q[1];
sx q[1];
rz(-2.5963929) q[1];
sx q[1];
rz(-1.6166519) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.10337457) q[0];
sx q[0];
rz(-0.056110121) q[0];
sx q[0];
rz(2.9655064) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.97648804) q[2];
sx q[2];
rz(-0.39820489) q[2];
sx q[2];
rz(1.9081685) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.9774919) q[1];
sx q[1];
rz(-2.3445446) q[1];
sx q[1];
rz(-2.1063095) q[1];
rz(-pi) q[2];
x q[2];
rz(2.3297263) q[3];
sx q[3];
rz(-0.97211876) q[3];
sx q[3];
rz(2.6177399) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-1.5357431) q[2];
sx q[2];
rz(-2.3652786) q[2];
sx q[2];
rz(-0.95174754) q[2];
rz(0.4176628) q[3];
sx q[3];
rz(-2.8916736) q[3];
sx q[3];
rz(1.1550268) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8010537) q[0];
sx q[0];
rz(-1.2507573) q[0];
sx q[0];
rz(-0.75403768) q[0];
rz(2.2484089) q[1];
sx q[1];
rz(-1.040753) q[1];
sx q[1];
rz(0.16597861) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6245859) q[0];
sx q[0];
rz(-0.69147325) q[0];
sx q[0];
rz(2.9777479) q[0];
rz(-pi) q[1];
x q[1];
rz(1.3909666) q[2];
sx q[2];
rz(-1.0713312) q[2];
sx q[2];
rz(-1.3701554) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-0.41478911) q[1];
sx q[1];
rz(-0.94894743) q[1];
sx q[1];
rz(-1.0110823) q[1];
x q[2];
rz(1.2955722) q[3];
sx q[3];
rz(-1.7050192) q[3];
sx q[3];
rz(-1.272066) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-0.37504998) q[2];
sx q[2];
rz(-2.0039717) q[2];
sx q[2];
rz(0.005391187) q[2];
rz(-3.052875) q[3];
sx q[3];
rz(-1.5327449) q[3];
sx q[3];
rz(0.79571342) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8832815) q[0];
sx q[0];
rz(-1.9323876) q[0];
sx q[0];
rz(-0.2051556) q[0];
rz(-0.908665) q[1];
sx q[1];
rz(-2.2712207) q[1];
sx q[1];
rz(-0.044513449) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9712898) q[0];
sx q[0];
rz(-1.4409541) q[0];
sx q[0];
rz(1.4128039) q[0];
rz(-pi) q[1];
rz(1.3088063) q[2];
sx q[2];
rz(-0.14523187) q[2];
sx q[2];
rz(-1.1373718) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(1.2878814) q[1];
sx q[1];
rz(-1.7856303) q[1];
sx q[1];
rz(-0.98835215) q[1];
rz(1.3839339) q[3];
sx q[3];
rz(-3.0198041) q[3];
sx q[3];
rz(0.78021061) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-2.8160416) q[2];
sx q[2];
rz(-2.6191235) q[2];
sx q[2];
rz(-2.6204056) q[2];
rz(2.9442287) q[3];
sx q[3];
rz(-1.3857931) q[3];
sx q[3];
rz(-0.50281966) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
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
rz(0.40081438) q[0];
sx q[0];
rz(-2.9476808) q[0];
sx q[0];
rz(0.25522301) q[0];
rz(0.1420282) q[1];
sx q[1];
rz(-1.4165001) q[1];
sx q[1];
rz(-0.35370383) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2784894) q[0];
sx q[0];
rz(-1.8971839) q[0];
sx q[0];
rz(1.2005689) q[0];
rz(-pi) q[1];
x q[1];
rz(0.50019294) q[2];
sx q[2];
rz(-1.0076773) q[2];
sx q[2];
rz(-1.4266071) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(0.6981909) q[1];
sx q[1];
rz(-2.576722) q[1];
sx q[1];
rz(-0.90157149) q[1];
x q[2];
rz(2.5442637) q[3];
sx q[3];
rz(-2.2976365) q[3];
sx q[3];
rz(-2.8387808) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.8346617) q[2];
sx q[2];
rz(-2.8262409) q[2];
sx q[2];
rz(-1.5241415) q[2];
rz(2.2751685) q[3];
sx q[3];
rz(-1.6774991) q[3];
sx q[3];
rz(-2.5518937) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.50188142) q[0];
sx q[0];
rz(-2.9852133) q[0];
sx q[0];
rz(1.7891275) q[0];
rz(1.9068708) q[1];
sx q[1];
rz(-2.9094978) q[1];
sx q[1];
rz(-2.629705) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5049625) q[0];
sx q[0];
rz(-1.413336) q[0];
sx q[0];
rz(0.053863346) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.9779766) q[2];
sx q[2];
rz(-2.3694042) q[2];
sx q[2];
rz(-0.26547394) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-1.3583957) q[1];
sx q[1];
rz(-1.277521) q[1];
sx q[1];
rz(2.6986928) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.98661042) q[3];
sx q[3];
rz(-1.8711149) q[3];
sx q[3];
rz(-2.9943313) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-1.1341683) q[2];
sx q[2];
rz(-1.3205426) q[2];
sx q[2];
rz(2.7443938) q[2];
rz(1.0154137) q[3];
sx q[3];
rz(-2.1404603) q[3];
sx q[3];
rz(2.5889682) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7427202) q[0];
sx q[0];
rz(-1.9341368) q[0];
sx q[0];
rz(0.43750986) q[0];
rz(-0.73879009) q[1];
sx q[1];
rz(-2.3414325) q[1];
sx q[1];
rz(-1.4791666) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2302983) q[0];
sx q[0];
rz(-1.084096) q[0];
sx q[0];
rz(-1.507797) q[0];
rz(-pi) q[1];
x q[1];
rz(3.0257341) q[2];
sx q[2];
rz(-1.9278952) q[2];
sx q[2];
rz(3.0810205) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(0.38592713) q[1];
sx q[1];
rz(-0.12317056) q[1];
sx q[1];
rz(2.6120899) q[1];
rz(0.86410256) q[3];
sx q[3];
rz(-2.2165944) q[3];
sx q[3];
rz(-1.0458664) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.9208357) q[2];
sx q[2];
rz(-1.3584542) q[2];
sx q[2];
rz(-0.1753359) q[2];
rz(-2.1389424) q[3];
sx q[3];
rz(-0.23313871) q[3];
sx q[3];
rz(-1.233915) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.47068448) q[0];
sx q[0];
rz(-1.4574454) q[0];
sx q[0];
rz(-1.1048143) q[0];
rz(0.15057527) q[1];
sx q[1];
rz(-0.87599788) q[1];
sx q[1];
rz(-1.8465975) q[1];
rz(2.9433123) q[2];
sx q[2];
rz(-1.5783327) q[2];
sx q[2];
rz(-0.21680149) q[2];
rz(-0.2507052) q[3];
sx q[3];
rz(-1.7964994) q[3];
sx q[3];
rz(-2.9356706) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
