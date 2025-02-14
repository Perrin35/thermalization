OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(0.72260296) q[0];
sx q[0];
rz(-2.498772) q[0];
sx q[0];
rz(8.2052054) q[0];
rz(0.83528432) q[1];
sx q[1];
rz(-1.3568027) q[1];
sx q[1];
rz(-2.2396483) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.038739) q[0];
sx q[0];
rz(-1.7243288) q[0];
sx q[0];
rz(0.6153591) q[0];
rz(-1.7576394) q[2];
sx q[2];
rz(-1.1964436) q[2];
sx q[2];
rz(-0.49546212) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(0.28022596) q[1];
sx q[1];
rz(-1.8609443) q[1];
sx q[1];
rz(2.9600083) q[1];
rz(-pi) q[2];
rz(0.70848453) q[3];
sx q[3];
rz(-0.86743858) q[3];
sx q[3];
rz(0.10242232) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-2.1787662) q[2];
sx q[2];
rz(-1.0054532) q[2];
sx q[2];
rz(-0.82628769) q[2];
rz(-1.7360342) q[3];
sx q[3];
rz(-2.8345351) q[3];
sx q[3];
rz(-0.74339408) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.58251441) q[0];
sx q[0];
rz(-2.7358416) q[0];
sx q[0];
rz(-1.7412809) q[0];
rz(2.0179613) q[1];
sx q[1];
rz(-2.7476937) q[1];
sx q[1];
rz(-1.0981015) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0880249) q[0];
sx q[0];
rz(-2.9505749) q[0];
sx q[0];
rz(2.5294526) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.2400572) q[2];
sx q[2];
rz(-1.8914701) q[2];
sx q[2];
rz(1.641524) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(0.73126331) q[1];
sx q[1];
rz(-1.6917233) q[1];
sx q[1];
rz(1.4951453) q[1];
rz(-0.16657515) q[3];
sx q[3];
rz(-1.9431356) q[3];
sx q[3];
rz(0.53322116) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(2.0565679) q[2];
sx q[2];
rz(-2.8242064) q[2];
sx q[2];
rz(1.8918096) q[2];
rz(-1.362484) q[3];
sx q[3];
rz(-1.0008078) q[3];
sx q[3];
rz(-1.6574297) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
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
rz(-0.63610858) q[0];
sx q[0];
rz(-0.69378575) q[0];
sx q[0];
rz(-2.8662477) q[0];
rz(2.6823726) q[1];
sx q[1];
rz(-0.95454916) q[1];
sx q[1];
rz(-0.053344639) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.1090272) q[0];
sx q[0];
rz(-1.7836187) q[0];
sx q[0];
rz(-1.4906917) q[0];
rz(-3.1231327) q[2];
sx q[2];
rz(-2.0316796) q[2];
sx q[2];
rz(0.26319749) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-0.054823067) q[1];
sx q[1];
rz(-2.5882698) q[1];
sx q[1];
rz(0.077484681) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.8998013) q[3];
sx q[3];
rz(-1.3734372) q[3];
sx q[3];
rz(-0.98361919) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-0.82103819) q[2];
sx q[2];
rz(-2.7802763) q[2];
sx q[2];
rz(-1.6374755) q[2];
rz(1.5004246) q[3];
sx q[3];
rz(-1.0502366) q[3];
sx q[3];
rz(-2.8659081) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.42030537) q[0];
sx q[0];
rz(-0.5834226) q[0];
sx q[0];
rz(2.6570008) q[0];
rz(-1.8991607) q[1];
sx q[1];
rz(-2.395605) q[1];
sx q[1];
rz(0.34908435) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9681681) q[0];
sx q[0];
rz(-1.613488) q[0];
sx q[0];
rz(2.1195565) q[0];
rz(-pi) q[1];
rz(-1.3626839) q[2];
sx q[2];
rz(-0.35018626) q[2];
sx q[2];
rz(-1.1454358) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-3.1368833) q[1];
sx q[1];
rz(-1.3420958) q[1];
sx q[1];
rz(1.7654224) q[1];
rz(-pi) q[2];
x q[2];
rz(2.1717668) q[3];
sx q[3];
rz(-1.2194014) q[3];
sx q[3];
rz(-1.9072208) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(2.7829973) q[2];
sx q[2];
rz(-1.2025669) q[2];
sx q[2];
rz(-0.45004582) q[2];
rz(2.3027244) q[3];
sx q[3];
rz(-1.5476371) q[3];
sx q[3];
rz(-0.1990327) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.49804509) q[0];
sx q[0];
rz(-1.4692551) q[0];
sx q[0];
rz(3.041748) q[0];
rz(-1.8210583) q[1];
sx q[1];
rz(-2.1456199) q[1];
sx q[1];
rz(-2.5340396) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.71220416) q[0];
sx q[0];
rz(-0.11183248) q[0];
sx q[0];
rz(2.4080316) q[0];
rz(-pi) q[1];
rz(2.8855912) q[2];
sx q[2];
rz(-0.56016541) q[2];
sx q[2];
rz(-1.7988811) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-2.1888426) q[1];
sx q[1];
rz(-2.5217068) q[1];
sx q[1];
rz(-0.37921885) q[1];
rz(1.3123061) q[3];
sx q[3];
rz(-1.6655169) q[3];
sx q[3];
rz(-2.3324089) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-1.1943835) q[2];
sx q[2];
rz(-0.76346976) q[2];
sx q[2];
rz(-2.6109931) q[2];
rz(-2.7001906) q[3];
sx q[3];
rz(-2.1680021) q[3];
sx q[3];
rz(2.1921564) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.43668231) q[0];
sx q[0];
rz(-1.7969776) q[0];
sx q[0];
rz(0.5740903) q[0];
rz(-0.87999815) q[1];
sx q[1];
rz(-2.0206385) q[1];
sx q[1];
rz(-1.3425739) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8594151) q[0];
sx q[0];
rz(-1.1033774) q[0];
sx q[0];
rz(-0.22478454) q[0];
x q[1];
rz(0.75127496) q[2];
sx q[2];
rz(-0.57005586) q[2];
sx q[2];
rz(1.3700203) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-1.3664535) q[1];
sx q[1];
rz(-0.43700686) q[1];
sx q[1];
rz(-0.98458146) q[1];
x q[2];
rz(2.4785512) q[3];
sx q[3];
rz(-1.6812857) q[3];
sx q[3];
rz(-0.35943951) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(2.4673956) q[2];
sx q[2];
rz(-0.28417045) q[2];
sx q[2];
rz(-2.6944842) q[2];
rz(-1.0304662) q[3];
sx q[3];
rz(-1.5117398) q[3];
sx q[3];
rz(-2.7569421) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.36987385) q[0];
sx q[0];
rz(-1.2971224) q[0];
sx q[0];
rz(2.7287927) q[0];
rz(-2.0665456) q[1];
sx q[1];
rz(-2.0037035) q[1];
sx q[1];
rz(-0.80361754) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5301054) q[0];
sx q[0];
rz(-0.75615935) q[0];
sx q[0];
rz(2.3458781) q[0];
rz(-pi) q[1];
rz(-2.2718524) q[2];
sx q[2];
rz(-0.77574965) q[2];
sx q[2];
rz(2.9091331) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-0.072783006) q[1];
sx q[1];
rz(-2.3549665) q[1];
sx q[1];
rz(-1.5907611) q[1];
rz(-pi) q[2];
rz(-0.39866205) q[3];
sx q[3];
rz(-2.9473262) q[3];
sx q[3];
rz(0.27856058) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-0.86218086) q[2];
sx q[2];
rz(-1.923424) q[2];
sx q[2];
rz(1.5009521) q[2];
rz(2.5189279) q[3];
sx q[3];
rz(-0.77087918) q[3];
sx q[3];
rz(-1.6525432) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.415446) q[0];
sx q[0];
rz(-1.5203238) q[0];
sx q[0];
rz(2.5836482) q[0];
rz(-0.56124148) q[1];
sx q[1];
rz(-1.3713501) q[1];
sx q[1];
rz(-1.7254613) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5197882) q[0];
sx q[0];
rz(-2.0178231) q[0];
sx q[0];
rz(0.31479737) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.1886732) q[2];
sx q[2];
rz(-1.2968569) q[2];
sx q[2];
rz(0.94715202) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.76991612) q[1];
sx q[1];
rz(-1.8990108) q[1];
sx q[1];
rz(-0.10336831) q[1];
rz(-pi) q[2];
rz(-0.21702311) q[3];
sx q[3];
rz(-0.76730928) q[3];
sx q[3];
rz(1.2485152) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(2.2531565) q[2];
sx q[2];
rz(-2.3641868) q[2];
sx q[2];
rz(-2.1412795) q[2];
rz(0.12217626) q[3];
sx q[3];
rz(-1.6413942) q[3];
sx q[3];
rz(0.15672556) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.67824739) q[0];
sx q[0];
rz(-1.1868917) q[0];
sx q[0];
rz(3.1372702) q[0];
rz(-0.16383544) q[1];
sx q[1];
rz(-2.7163353) q[1];
sx q[1];
rz(1.3660627) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3177494) q[0];
sx q[0];
rz(-1.5986406) q[0];
sx q[0];
rz(-2.2180024) q[0];
rz(-pi) q[1];
rz(1.2605569) q[2];
sx q[2];
rz(-1.2149335) q[2];
sx q[2];
rz(-2.9663393) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-2.8722265) q[1];
sx q[1];
rz(-2.8878643) q[1];
sx q[1];
rz(0.14464186) q[1];
rz(-2.4716464) q[3];
sx q[3];
rz(-2.1996227) q[3];
sx q[3];
rz(-0.82312246) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-2.0020478) q[2];
sx q[2];
rz(-0.9674955) q[2];
sx q[2];
rz(0.61919332) q[2];
rz(0.55245095) q[3];
sx q[3];
rz(-1.8360454) q[3];
sx q[3];
rz(-2.9964871) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7904952) q[0];
sx q[0];
rz(-2.9572697) q[0];
sx q[0];
rz(0.47875324) q[0];
rz(1.152773) q[1];
sx q[1];
rz(-0.29007998) q[1];
sx q[1];
rz(-0.94720381) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9298332) q[0];
sx q[0];
rz(-2.9538028) q[0];
sx q[0];
rz(-1.7540356) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.2113153) q[2];
sx q[2];
rz(-0.61865846) q[2];
sx q[2];
rz(1.930069) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-2.3403907) q[1];
sx q[1];
rz(-1.8762734) q[1];
sx q[1];
rz(0.82251541) q[1];
rz(-pi) q[2];
rz(0.35263108) q[3];
sx q[3];
rz(-1.8049587) q[3];
sx q[3];
rz(1.4591097) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(1.3843711) q[2];
sx q[2];
rz(-0.20558509) q[2];
sx q[2];
rz(-2.8749386) q[2];
rz(1.0673149) q[3];
sx q[3];
rz(-1.9284733) q[3];
sx q[3];
rz(-0.96323577) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
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
rz(-0.08854475) q[0];
sx q[0];
rz(-1.5374669) q[0];
sx q[0];
rz(1.5446825) q[0];
rz(2.0883941) q[1];
sx q[1];
rz(-1.2146626) q[1];
sx q[1];
rz(0.76520898) q[1];
rz(-3.13022) q[2];
sx q[2];
rz(-0.24200242) q[2];
sx q[2];
rz(0.92462362) q[2];
rz(-2.8070634) q[3];
sx q[3];
rz(-2.5361037) q[3];
sx q[3];
rz(-1.1246214) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
