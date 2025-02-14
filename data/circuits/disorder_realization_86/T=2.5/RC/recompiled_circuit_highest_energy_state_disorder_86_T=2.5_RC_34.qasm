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
rz(0.64492172) q[0];
sx q[0];
rz(-0.46907297) q[0];
sx q[0];
rz(2.1512349) q[0];
rz(-1.1733836) q[1];
sx q[1];
rz(2.2819509) q[1];
sx q[1];
rz(11.717164) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9678803) q[0];
sx q[0];
rz(-0.033941887) q[0];
sx q[0];
rz(1.6544266) q[0];
x q[1];
rz(-2.8013632) q[2];
sx q[2];
rz(-2.2417667) q[2];
sx q[2];
rz(-3.1052507) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-1.7265004) q[1];
sx q[1];
rz(-2.1884568) q[1];
sx q[1];
rz(0.32869494) q[1];
rz(1.0484656) q[3];
sx q[3];
rz(-1.4712787) q[3];
sx q[3];
rz(-2.8985648) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-2.2363756) q[2];
sx q[2];
rz(-1.507501) q[2];
sx q[2];
rz(-0.53984731) q[2];
rz(-1.5652462) q[3];
sx q[3];
rz(-2.7326475) q[3];
sx q[3];
rz(-1.7319771) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0419615) q[0];
sx q[0];
rz(-2.4653682) q[0];
sx q[0];
rz(1.3270295) q[0];
rz(2.3565893) q[1];
sx q[1];
rz(-1.4229341) q[1];
sx q[1];
rz(0.50055093) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6645522) q[0];
sx q[0];
rz(-1.0196166) q[0];
sx q[0];
rz(-1.5413947) q[0];
rz(-pi) q[1];
x q[1];
rz(2.9517216) q[2];
sx q[2];
rz(-1.9315757) q[2];
sx q[2];
rz(-1.2072414) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(3.0620898) q[1];
sx q[1];
rz(-2.4356053) q[1];
sx q[1];
rz(2.1022309) q[1];
x q[2];
rz(-0.80970069) q[3];
sx q[3];
rz(-1.9040344) q[3];
sx q[3];
rz(-2.7908195) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(2.0565722) q[2];
sx q[2];
rz(-2.4041921) q[2];
sx q[2];
rz(2.6596587) q[2];
rz(2.7815172) q[3];
sx q[3];
rz(-1.998338) q[3];
sx q[3];
rz(2.9329407) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8423186) q[0];
sx q[0];
rz(-2.0330918) q[0];
sx q[0];
rz(-0.47766787) q[0];
rz(2.281669) q[1];
sx q[1];
rz(-1.0483024) q[1];
sx q[1];
rz(0.28712505) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5764113) q[0];
sx q[0];
rz(-0.14248304) q[0];
sx q[0];
rz(-1.0866685) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.3541917) q[2];
sx q[2];
rz(-1.9176716) q[2];
sx q[2];
rz(1.5225449) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.9707253) q[1];
sx q[1];
rz(-1.5795465) q[1];
sx q[1];
rz(0.88717242) q[1];
rz(-pi) q[2];
x q[2];
rz(-3.0564553) q[3];
sx q[3];
rz(-1.9508024) q[3];
sx q[3];
rz(0.054084965) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-2.3890248) q[2];
sx q[2];
rz(-1.7914881) q[2];
sx q[2];
rz(-3.1217421) q[2];
rz(-2.1974473) q[3];
sx q[3];
rz(-2.4343334) q[3];
sx q[3];
rz(2.7437362) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.95553628) q[0];
sx q[0];
rz(-0.63025403) q[0];
sx q[0];
rz(1.3007042) q[0];
rz(-2.6553254) q[1];
sx q[1];
rz(-1.9673037) q[1];
sx q[1];
rz(3.1062612) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6636642) q[0];
sx q[0];
rz(-1.999117) q[0];
sx q[0];
rz(0.028686319) q[0];
rz(-2.7443419) q[2];
sx q[2];
rz(-1.78671) q[2];
sx q[2];
rz(-1.6524894) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-1.1580646) q[1];
sx q[1];
rz(-3.1175545) q[1];
sx q[1];
rz(2.9062494) q[1];
rz(-0.70628665) q[3];
sx q[3];
rz(-2.379198) q[3];
sx q[3];
rz(-1.958999) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(2.6639634) q[2];
sx q[2];
rz(-0.56325459) q[2];
sx q[2];
rz(0.25838724) q[2];
rz(-0.7102617) q[3];
sx q[3];
rz(-0.60451549) q[3];
sx q[3];
rz(1.5593504) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.48359394) q[0];
sx q[0];
rz(-2.3778264) q[0];
sx q[0];
rz(-2.4972231) q[0];
rz(-0.98980347) q[1];
sx q[1];
rz(-0.80392307) q[1];
sx q[1];
rz(2.03233) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8744226) q[0];
sx q[0];
rz(-2.5760898) q[0];
sx q[0];
rz(1.1995951) q[0];
rz(2.0844056) q[2];
sx q[2];
rz(-1.6053146) q[2];
sx q[2];
rz(-2.2396954) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(0.39711827) q[1];
sx q[1];
rz(-2.3496685) q[1];
sx q[1];
rz(0.20361118) q[1];
rz(1.8568138) q[3];
sx q[3];
rz(-0.83854616) q[3];
sx q[3];
rz(0.58457182) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-0.41205078) q[2];
sx q[2];
rz(-1.5645626) q[2];
sx q[2];
rz(-2.5315419) q[2];
rz(3.0610906) q[3];
sx q[3];
rz(-0.1463612) q[3];
sx q[3];
rz(3.1350873) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.42008156) q[0];
sx q[0];
rz(-0.93669909) q[0];
sx q[0];
rz(-1.6625846) q[0];
rz(-1.1889907) q[1];
sx q[1];
rz(-2.7956796) q[1];
sx q[1];
rz(1.4422653) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8540031) q[0];
sx q[0];
rz(-0.95335273) q[0];
sx q[0];
rz(-2.94728) q[0];
rz(-pi) q[1];
rz(-2.5601588) q[2];
sx q[2];
rz(-1.5155041) q[2];
sx q[2];
rz(-0.96382574) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(0.4870781) q[1];
sx q[1];
rz(-2.2527825) q[1];
sx q[1];
rz(2.712516) q[1];
x q[2];
rz(-0.54338065) q[3];
sx q[3];
rz(-0.60688775) q[3];
sx q[3];
rz(0.010963765) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-0.46001616) q[2];
sx q[2];
rz(-0.67395335) q[2];
sx q[2];
rz(-2.4665534) q[2];
rz(0.33440822) q[3];
sx q[3];
rz(-1.6655917) q[3];
sx q[3];
rz(0.36791754) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6595031) q[0];
sx q[0];
rz(-2.1559494) q[0];
sx q[0];
rz(3.0159045) q[0];
rz(-2.9478759) q[1];
sx q[1];
rz(-2.2592762) q[1];
sx q[1];
rz(1.9901336) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.73456681) q[0];
sx q[0];
rz(-2.1687963) q[0];
sx q[0];
rz(0.96877386) q[0];
rz(-pi) q[1];
x q[1];
rz(2.2073484) q[2];
sx q[2];
rz(-2.6912874) q[2];
sx q[2];
rz(0.59666419) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-2.0380504) q[1];
sx q[1];
rz(-2.0972431) q[1];
sx q[1];
rz(-1.0564694) q[1];
rz(-2.4667593) q[3];
sx q[3];
rz(-1.1911567) q[3];
sx q[3];
rz(1.4745431) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-1.2842399) q[2];
sx q[2];
rz(-1.7809296) q[2];
sx q[2];
rz(2.8947158) q[2];
rz(0.85865584) q[3];
sx q[3];
rz(-2.1571428) q[3];
sx q[3];
rz(0.27913678) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.7849279) q[0];
sx q[0];
rz(-0.44342884) q[0];
sx q[0];
rz(0.32522935) q[0];
rz(-0.52109703) q[1];
sx q[1];
rz(-0.63322133) q[1];
sx q[1];
rz(-0.31347832) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3081261) q[0];
sx q[0];
rz(-2.0692056) q[0];
sx q[0];
rz(0.68584401) q[0];
rz(-pi) q[1];
rz(-3.0738906) q[2];
sx q[2];
rz(-1.3460068) q[2];
sx q[2];
rz(2.1613286) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-1.5969475) q[1];
sx q[1];
rz(-1.2352422) q[1];
sx q[1];
rz(0.38864522) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.9133592) q[3];
sx q[3];
rz(-0.54751626) q[3];
sx q[3];
rz(-1.3207796) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(0.43361214) q[2];
sx q[2];
rz(-1.2182451) q[2];
sx q[2];
rz(1.8617967) q[2];
rz(-0.72569877) q[3];
sx q[3];
rz(-1.9819219) q[3];
sx q[3];
rz(-2.3316135) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5028266) q[0];
sx q[0];
rz(-1.1932729) q[0];
sx q[0];
rz(2.9916812) q[0];
rz(1.8984849) q[1];
sx q[1];
rz(-1.5877692) q[1];
sx q[1];
rz(1.6601723) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7090764) q[0];
sx q[0];
rz(-3.0750599) q[0];
sx q[0];
rz(-0.82609047) q[0];
rz(-pi) q[1];
x q[1];
rz(0.051542087) q[2];
sx q[2];
rz(-1.3304425) q[2];
sx q[2];
rz(0.86567825) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.77566389) q[1];
sx q[1];
rz(-1.3559113) q[1];
sx q[1];
rz(0.38916969) q[1];
x q[2];
rz(1.0970988) q[3];
sx q[3];
rz(-0.27333958) q[3];
sx q[3];
rz(-2.7784082) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-2.3159065) q[2];
sx q[2];
rz(-1.76182) q[2];
sx q[2];
rz(-2.1659577) q[2];
rz(-2.0129096) q[3];
sx q[3];
rz(-2.5895139) q[3];
sx q[3];
rz(2.8868207) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.16633701) q[0];
sx q[0];
rz(-1.0003426) q[0];
sx q[0];
rz(0.16217232) q[0];
rz(1.8099248) q[1];
sx q[1];
rz(-2.5089896) q[1];
sx q[1];
rz(0.046028927) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8019288) q[0];
sx q[0];
rz(-1.5886512) q[0];
sx q[0];
rz(0.0059203832) q[0];
x q[1];
rz(0.46398766) q[2];
sx q[2];
rz(-1.2210238) q[2];
sx q[2];
rz(-0.9198364) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(1.6190336) q[1];
sx q[1];
rz(-1.7184034) q[1];
sx q[1];
rz(-1.1655432) q[1];
rz(-1.3999942) q[3];
sx q[3];
rz(-1.9055589) q[3];
sx q[3];
rz(-0.46250175) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(0.33596805) q[2];
sx q[2];
rz(-0.38463548) q[2];
sx q[2];
rz(-1.1495205) q[2];
rz(2.6286821) q[3];
sx q[3];
rz(-1.3866813) q[3];
sx q[3];
rz(-2.8368867) q[3];
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
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.49160663) q[0];
sx q[0];
rz(-2.2254324) q[0];
sx q[0];
rz(-1.9944763) q[0];
rz(0.06123771) q[1];
sx q[1];
rz(-1.9495268) q[1];
sx q[1];
rz(-1.3757642) q[1];
rz(-2.1381151) q[2];
sx q[2];
rz(-0.85366953) q[2];
sx q[2];
rz(0.015346957) q[2];
rz(0.41027222) q[3];
sx q[3];
rz(-0.82220746) q[3];
sx q[3];
rz(-2.0668277) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
