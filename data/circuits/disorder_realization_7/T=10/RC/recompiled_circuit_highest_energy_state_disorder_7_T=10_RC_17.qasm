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
rz(1.3697019) q[0];
sx q[0];
rz(-0.43039027) q[0];
sx q[0];
rz(-0.16464591) q[0];
rz(-1.6904263) q[1];
sx q[1];
rz(-0.37835205) q[1];
sx q[1];
rz(-2.4317256) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9116316) q[0];
sx q[0];
rz(-0.29567406) q[0];
sx q[0];
rz(-1.1082746) q[0];
x q[1];
rz(0.69605445) q[2];
sx q[2];
rz(-0.6185607) q[2];
sx q[2];
rz(1.6060917) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(1.134608) q[1];
sx q[1];
rz(-0.35672327) q[1];
sx q[1];
rz(-1.7118042) q[1];
rz(0.23386896) q[3];
sx q[3];
rz(-1.4483671) q[3];
sx q[3];
rz(2.8120638) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-1.5300488) q[2];
sx q[2];
rz(-2.7908466) q[2];
sx q[2];
rz(1.9225527) q[2];
rz(-0.48398316) q[3];
sx q[3];
rz(-2.2958906) q[3];
sx q[3];
rz(1.774196) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.71500635) q[0];
sx q[0];
rz(-0.33892092) q[0];
sx q[0];
rz(1.4434848) q[0];
rz(-1.8679856) q[1];
sx q[1];
rz(-2.1111635) q[1];
sx q[1];
rz(0.66676203) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6359207) q[0];
sx q[0];
rz(-1.8673602) q[0];
sx q[0];
rz(-1.7478329) q[0];
rz(-pi) q[1];
rz(-0.86033852) q[2];
sx q[2];
rz(-0.14709148) q[2];
sx q[2];
rz(1.8837613) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-2.8005008) q[1];
sx q[1];
rz(-1.9750764) q[1];
sx q[1];
rz(2.0501818) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.1902892) q[3];
sx q[3];
rz(-2.3535121) q[3];
sx q[3];
rz(0.57315592) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-2.4296809) q[2];
sx q[2];
rz(-1.9671755) q[2];
sx q[2];
rz(-2.5453117) q[2];
rz(1.0239673) q[3];
sx q[3];
rz(-1.3919316) q[3];
sx q[3];
rz(-2.021324) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4100818) q[0];
sx q[0];
rz(-0.54786587) q[0];
sx q[0];
rz(-1.6743073) q[0];
rz(0.038854988) q[1];
sx q[1];
rz(-1.7170186) q[1];
sx q[1];
rz(-0.88599667) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.17490921) q[0];
sx q[0];
rz(-1.0499448) q[0];
sx q[0];
rz(1.0731927) q[0];
x q[1];
rz(1.5019234) q[2];
sx q[2];
rz(-0.67126545) q[2];
sx q[2];
rz(-1.4466937) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-1.6941144) q[1];
sx q[1];
rz(-1.9285818) q[1];
sx q[1];
rz(1.5286137) q[1];
x q[2];
rz(-1.6502077) q[3];
sx q[3];
rz(-1.5089499) q[3];
sx q[3];
rz(-1.1140149) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(0.92750612) q[2];
sx q[2];
rz(-1.4518041) q[2];
sx q[2];
rz(-0.56373325) q[2];
rz(0.8756513) q[3];
sx q[3];
rz(-1.0522269) q[3];
sx q[3];
rz(-1.0912857) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.57986528) q[0];
sx q[0];
rz(-2.2932597) q[0];
sx q[0];
rz(2.0735829) q[0];
rz(-0.39455286) q[1];
sx q[1];
rz(-2.6119472) q[1];
sx q[1];
rz(1.4716757) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.040634) q[0];
sx q[0];
rz(-1.9143595) q[0];
sx q[0];
rz(-0.45303194) q[0];
rz(-pi) q[1];
x q[1];
rz(1.8561697) q[2];
sx q[2];
rz(-1.4625408) q[2];
sx q[2];
rz(-1.4982323) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(0.79395548) q[1];
sx q[1];
rz(-0.82073602) q[1];
sx q[1];
rz(2.1669037) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.25217523) q[3];
sx q[3];
rz(-2.7311595) q[3];
sx q[3];
rz(1.3595734) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(2.7431405) q[2];
sx q[2];
rz(-1.8299711) q[2];
sx q[2];
rz(3.0121646) q[2];
rz(-0.5433003) q[3];
sx q[3];
rz(-1.0933417) q[3];
sx q[3];
rz(1.7530542) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1035136) q[0];
sx q[0];
rz(-1.0032126) q[0];
sx q[0];
rz(0.45573086) q[0];
rz(0.56574667) q[1];
sx q[1];
rz(-2.5028298) q[1];
sx q[1];
rz(-2.9260213) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5181464) q[0];
sx q[0];
rz(-1.7404273) q[0];
sx q[0];
rz(0.87624936) q[0];
rz(-pi) q[1];
x q[1];
rz(2.2224765) q[2];
sx q[2];
rz(-2.060215) q[2];
sx q[2];
rz(2.8280743) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.6681155) q[1];
sx q[1];
rz(-1.6398506) q[1];
sx q[1];
rz(-1.3051239) q[1];
x q[2];
rz(-3.0318063) q[3];
sx q[3];
rz(-1.2874787) q[3];
sx q[3];
rz(-2.7840419) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-2.4244708) q[2];
sx q[2];
rz(-2.8346546) q[2];
sx q[2];
rz(1.7090939) q[2];
rz(-2.0004382) q[3];
sx q[3];
rz(-1.9211831) q[3];
sx q[3];
rz(2.6728163) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5070709) q[0];
sx q[0];
rz(-2.2880726) q[0];
sx q[0];
rz(-2.5216907) q[0];
rz(-1.2076521) q[1];
sx q[1];
rz(-2.4687605) q[1];
sx q[1];
rz(2.8501453) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0196592) q[0];
sx q[0];
rz(-1.013151) q[0];
sx q[0];
rz(0.19361051) q[0];
rz(-pi) q[1];
rz(2.8180505) q[2];
sx q[2];
rz(-1.4446044) q[2];
sx q[2];
rz(3.0033811) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.92293948) q[1];
sx q[1];
rz(-2.1497552) q[1];
sx q[1];
rz(1.8576937) q[1];
rz(1.5602371) q[3];
sx q[3];
rz(-2.227475) q[3];
sx q[3];
rz(2.9845723) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-2.7181299) q[2];
sx q[2];
rz(-2.4726157) q[2];
sx q[2];
rz(-2.1675229) q[2];
rz(1.2232716) q[3];
sx q[3];
rz(-1.7134824) q[3];
sx q[3];
rz(2.1035002) q[3];
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
sx q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.85422) q[0];
sx q[0];
rz(-1.1690451) q[0];
sx q[0];
rz(-0.30782345) q[0];
rz(-2.8616915) q[1];
sx q[1];
rz(-2.8925536) q[1];
sx q[1];
rz(-1.8797967) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.38378206) q[0];
sx q[0];
rz(-1.6047262) q[0];
sx q[0];
rz(-1.4365804) q[0];
rz(-pi) q[1];
x q[1];
rz(0.73195157) q[2];
sx q[2];
rz(-1.2748655) q[2];
sx q[2];
rz(1.0481121) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(1.4018462) q[1];
sx q[1];
rz(-1.8658499) q[1];
sx q[1];
rz(3.1298545) q[1];
x q[2];
rz(-0.070075017) q[3];
sx q[3];
rz(-2.0349906) q[3];
sx q[3];
rz(-1.3843368) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(1.5032924) q[2];
sx q[2];
rz(-1.952848) q[2];
sx q[2];
rz(-2.3390181) q[2];
rz(-1.59683) q[3];
sx q[3];
rz(-0.39512008) q[3];
sx q[3];
rz(-1.4748658) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.06642) q[0];
sx q[0];
rz(-1.5936699) q[0];
sx q[0];
rz(-2.5286034) q[0];
rz(-1.0523484) q[1];
sx q[1];
rz(-2.3520825) q[1];
sx q[1];
rz(1.8863511) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.3259473) q[0];
sx q[0];
rz(-1.1425352) q[0];
sx q[0];
rz(1.0897348) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.9238812) q[2];
sx q[2];
rz(-1.048511) q[2];
sx q[2];
rz(-0.051607121) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-0.44574983) q[1];
sx q[1];
rz(-0.90409213) q[1];
sx q[1];
rz(-0.95335828) q[1];
rz(-1.4358712) q[3];
sx q[3];
rz(-1.8308911) q[3];
sx q[3];
rz(1.2430199) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-2.5455948) q[2];
sx q[2];
rz(-2.2277446) q[2];
sx q[2];
rz(0.14860281) q[2];
rz(2.4871067) q[3];
sx q[3];
rz(-1.9446257) q[3];
sx q[3];
rz(1.436208) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[3];
rz(pi/2) q[3];
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
rz(0.29110903) q[0];
sx q[0];
rz(-1.608526) q[0];
sx q[0];
rz(-3.1403551) q[0];
rz(1.8908063) q[1];
sx q[1];
rz(-0.93235278) q[1];
sx q[1];
rz(-2.1900182) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5244999) q[0];
sx q[0];
rz(-1.6902335) q[0];
sx q[0];
rz(-2.9725084) q[0];
x q[1];
rz(-0.96828332) q[2];
sx q[2];
rz(-0.70455019) q[2];
sx q[2];
rz(0.3144484) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(2.7504362) q[1];
sx q[1];
rz(-1.5215148) q[1];
sx q[1];
rz(1.5261914) q[1];
x q[2];
rz(2.9195027) q[3];
sx q[3];
rz(-0.40910334) q[3];
sx q[3];
rz(-2.3698185) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.1452267) q[2];
sx q[2];
rz(-1.1684912) q[2];
sx q[2];
rz(1.036198) q[2];
rz(-2.2717617) q[3];
sx q[3];
rz(-1.4709604) q[3];
sx q[3];
rz(-0.74365348) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
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
rz(-pi) q[0];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0319801) q[0];
sx q[0];
rz(-2.5514422) q[0];
sx q[0];
rz(1.1085283) q[0];
rz(-2.176586) q[1];
sx q[1];
rz(-1.7644278) q[1];
sx q[1];
rz(-0.56934294) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8764502) q[0];
sx q[0];
rz(-0.56395774) q[0];
sx q[0];
rz(-2.7258858) q[0];
rz(2.0351719) q[2];
sx q[2];
rz(-1.2341675) q[2];
sx q[2];
rz(-2.6576633) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-0.18817338) q[1];
sx q[1];
rz(-1.0480289) q[1];
sx q[1];
rz(-3.1320577) q[1];
rz(-0.17579349) q[3];
sx q[3];
rz(-2.3758278) q[3];
sx q[3];
rz(-0.24336963) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(1.0137332) q[2];
sx q[2];
rz(-0.49488417) q[2];
sx q[2];
rz(-0.046517046) q[2];
rz(0.30139309) q[3];
sx q[3];
rz(-1.2208341) q[3];
sx q[3];
rz(-0.2218328) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.90060577) q[0];
sx q[0];
rz(-1.2099246) q[0];
sx q[0];
rz(-1.8151617) q[0];
rz(-1.2251414) q[1];
sx q[1];
rz(-2.6381208) q[1];
sx q[1];
rz(-2.0874964) q[1];
rz(2.8711748) q[2];
sx q[2];
rz(-2.5106988) q[2];
sx q[2];
rz(-0.0028263447) q[2];
rz(-0.75792652) q[3];
sx q[3];
rz(-1.0316385) q[3];
sx q[3];
rz(1.3760174) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
