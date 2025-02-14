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
rz(1.6871356) q[0];
sx q[0];
rz(1.9782826) q[0];
sx q[0];
rz(10.692698) q[0];
rz(1.6549702) q[1];
sx q[1];
rz(-1.854719) q[1];
sx q[1];
rz(0.16965228) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9823224) q[0];
sx q[0];
rz(-1.9290315) q[0];
sx q[0];
rz(1.0591749) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.9965788) q[2];
sx q[2];
rz(-0.55934956) q[2];
sx q[2];
rz(-2.1970791) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(2.4625068) q[1];
sx q[1];
rz(-0.39606491) q[1];
sx q[1];
rz(-1.5157265) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.9291487) q[3];
sx q[3];
rz(-1.8768969) q[3];
sx q[3];
rz(-2.3116632) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(2.9156645) q[2];
sx q[2];
rz(-1.6253086) q[2];
sx q[2];
rz(0.020326745) q[2];
rz(2.9018371) q[3];
sx q[3];
rz(-2.8862947) q[3];
sx q[3];
rz(0.83893004) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
x q[2];
x q[3];
rz(pi/2) q[3];
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
rz(0.13424419) q[0];
sx q[0];
rz(-0.61779314) q[0];
sx q[0];
rz(-1.0963305) q[0];
rz(1.6657375) q[1];
sx q[1];
rz(-1.219039) q[1];
sx q[1];
rz(0.79442564) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6446895) q[0];
sx q[0];
rz(-2.6711426) q[0];
sx q[0];
rz(-0.75017922) q[0];
rz(0.75411039) q[2];
sx q[2];
rz(-2.0144785) q[2];
sx q[2];
rz(2.544683) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-2.8640722) q[1];
sx q[1];
rz(-1.9937859) q[1];
sx q[1];
rz(1.6865684) q[1];
rz(-pi) q[2];
rz(-1.7678472) q[3];
sx q[3];
rz(-2.3757977) q[3];
sx q[3];
rz(1.1546419) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(3.0558001) q[2];
sx q[2];
rz(-1.5689359) q[2];
sx q[2];
rz(-2.9158578) q[2];
rz(0.24383946) q[3];
sx q[3];
rz(-0.95832458) q[3];
sx q[3];
rz(2.9624511) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8826411) q[0];
sx q[0];
rz(-1.3132341) q[0];
sx q[0];
rz(-0.42695811) q[0];
rz(-0.34307617) q[1];
sx q[1];
rz(-1.9648569) q[1];
sx q[1];
rz(0.25161904) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.93676567) q[0];
sx q[0];
rz(-2.511095) q[0];
sx q[0];
rz(2.7555694) q[0];
x q[1];
rz(1.8015091) q[2];
sx q[2];
rz(-0.78015155) q[2];
sx q[2];
rz(3.0456269) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.15089825) q[1];
sx q[1];
rz(-2.7643235) q[1];
sx q[1];
rz(-2.6713085) q[1];
x q[2];
rz(1.512647) q[3];
sx q[3];
rz(-2.4302501) q[3];
sx q[3];
rz(2.1429246) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.130927) q[2];
sx q[2];
rz(-1.5074573) q[2];
sx q[2];
rz(-2.6962386) q[2];
rz(-1.917786) q[3];
sx q[3];
rz(-1.8755951) q[3];
sx q[3];
rz(0.38770097) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7447516) q[0];
sx q[0];
rz(-2.3362609) q[0];
sx q[0];
rz(-1.9824363) q[0];
rz(1.9388439) q[1];
sx q[1];
rz(-1.9614599) q[1];
sx q[1];
rz(2.4045827) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.979078) q[0];
sx q[0];
rz(-1.3231228) q[0];
sx q[0];
rz(-1.635627) q[0];
rz(-pi) q[1];
rz(2.7734408) q[2];
sx q[2];
rz(-0.92879229) q[2];
sx q[2];
rz(-2.8411691) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(1.678469) q[1];
sx q[1];
rz(-2.2003074) q[1];
sx q[1];
rz(1.786382) q[1];
x q[2];
rz(-2.8094711) q[3];
sx q[3];
rz(-0.17617036) q[3];
sx q[3];
rz(2.9834) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.1272588) q[2];
sx q[2];
rz(-0.88982439) q[2];
sx q[2];
rz(-2.3298402) q[2];
rz(-0.74770606) q[3];
sx q[3];
rz(-2.2816608) q[3];
sx q[3];
rz(0.060613304) q[3];
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
x q[0];
rz(-pi/2) q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.49374813) q[0];
sx q[0];
rz(-2.0525377) q[0];
sx q[0];
rz(-1.7895948) q[0];
rz(-0.79212517) q[1];
sx q[1];
rz(-2.242576) q[1];
sx q[1];
rz(2.1314714) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8787694) q[0];
sx q[0];
rz(-2.4973625) q[0];
sx q[0];
rz(-2.1488726) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.5158122) q[2];
sx q[2];
rz(-1.9628982) q[2];
sx q[2];
rz(0.34811172) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-1.5252851) q[1];
sx q[1];
rz(-2.1263564) q[1];
sx q[1];
rz(1.6643334) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.185818) q[3];
sx q[3];
rz(-2.9412651) q[3];
sx q[3];
rz(2.1708787) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(3.1048364) q[2];
sx q[2];
rz(-0.75883055) q[2];
sx q[2];
rz(-1.3834312) q[2];
rz(3.1116327) q[3];
sx q[3];
rz(-2.1432917) q[3];
sx q[3];
rz(-1.2768607) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6303915) q[0];
sx q[0];
rz(-0.41655219) q[0];
sx q[0];
rz(-0.1314441) q[0];
rz(0.99880544) q[1];
sx q[1];
rz(-1.9824332) q[1];
sx q[1];
rz(-1.3712032) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4192062) q[0];
sx q[0];
rz(-0.50438577) q[0];
sx q[0];
rz(-0.76865159) q[0];
x q[1];
rz(1.7741469) q[2];
sx q[2];
rz(-0.89471811) q[2];
sx q[2];
rz(-1.0695367) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(1.882521) q[1];
sx q[1];
rz(-1.4132199) q[1];
sx q[1];
rz(1.2941918) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.71089069) q[3];
sx q[3];
rz(-1.5734473) q[3];
sx q[3];
rz(-2.7972935) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(2.3580631) q[2];
sx q[2];
rz(-0.72856599) q[2];
sx q[2];
rz(1.0895458) q[2];
rz(0.36695668) q[3];
sx q[3];
rz(-0.4042545) q[3];
sx q[3];
rz(0.75016108) q[3];
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
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.628196) q[0];
sx q[0];
rz(-2.3371526) q[0];
sx q[0];
rz(-2.2688769) q[0];
rz(1.5645507) q[1];
sx q[1];
rz(-1.4955474) q[1];
sx q[1];
rz(-2.4094792) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.66595378) q[0];
sx q[0];
rz(-1.7717965) q[0];
sx q[0];
rz(2.4253286) q[0];
x q[1];
rz(-2.9896834) q[2];
sx q[2];
rz(-0.6912125) q[2];
sx q[2];
rz(0.37746261) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(2.8120576) q[1];
sx q[1];
rz(-2.0382216) q[1];
sx q[1];
rz(-0.17623623) q[1];
x q[2];
rz(1.7978396) q[3];
sx q[3];
rz(-1.4770835) q[3];
sx q[3];
rz(-2.7467487) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(2.6096036) q[2];
sx q[2];
rz(-2.1405818) q[2];
sx q[2];
rz(2.5606489) q[2];
rz(2.0161435) q[3];
sx q[3];
rz(-1.9476798) q[3];
sx q[3];
rz(-1.1552936) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4693212) q[0];
sx q[0];
rz(-3.0348365) q[0];
sx q[0];
rz(-2.971055) q[0];
rz(-1.2455617) q[1];
sx q[1];
rz(-1.6163369) q[1];
sx q[1];
rz(-0.68444288) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4218773) q[0];
sx q[0];
rz(-1.9027896) q[0];
sx q[0];
rz(-2.6396455) q[0];
rz(-pi) q[1];
rz(2.6587517) q[2];
sx q[2];
rz(-1.4582658) q[2];
sx q[2];
rz(-1.7140599) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(1.9065793) q[1];
sx q[1];
rz(-2.7827127) q[1];
sx q[1];
rz(-2.3770664) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.079250431) q[3];
sx q[3];
rz(-2.4911315) q[3];
sx q[3];
rz(-1.4440756) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(0.30274063) q[2];
sx q[2];
rz(-1.799823) q[2];
sx q[2];
rz(2.0212685) q[2];
rz(-0.66425792) q[3];
sx q[3];
rz(-2.7393326) q[3];
sx q[3];
rz(0.50814381) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6054194) q[0];
sx q[0];
rz(-2.7689731) q[0];
sx q[0];
rz(-2.5033409) q[0];
rz(3.1328746) q[1];
sx q[1];
rz(-1.1161085) q[1];
sx q[1];
rz(2.7353824) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2531949) q[0];
sx q[0];
rz(-2.8450845) q[0];
sx q[0];
rz(1.967708) q[0];
rz(-pi) q[1];
rz(1.137144) q[2];
sx q[2];
rz(-1.6182476) q[2];
sx q[2];
rz(-0.46772568) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(2.6510424) q[1];
sx q[1];
rz(-0.79572751) q[1];
sx q[1];
rz(2.9814238) q[1];
rz(-pi) q[2];
rz(1.5962192) q[3];
sx q[3];
rz(-1.4973728) q[3];
sx q[3];
rz(1.582445) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-0.093988769) q[2];
sx q[2];
rz(-1.1270019) q[2];
sx q[2];
rz(0.36063933) q[2];
rz(2.4141198) q[3];
sx q[3];
rz(-0.68220264) q[3];
sx q[3];
rz(2.8235161) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3280846) q[0];
sx q[0];
rz(-2.1961975) q[0];
sx q[0];
rz(2.9720921) q[0];
rz(0.70973474) q[1];
sx q[1];
rz(-1.7252445) q[1];
sx q[1];
rz(2.0555326) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5026131) q[0];
sx q[0];
rz(-1.7923681) q[0];
sx q[0];
rz(-2.6022439) q[0];
rz(-0.96699826) q[2];
sx q[2];
rz(-0.51568401) q[2];
sx q[2];
rz(-2.5623851) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.5063613) q[1];
sx q[1];
rz(-1.4725793) q[1];
sx q[1];
rz(2.1327095) q[1];
rz(-pi) q[2];
x q[2];
rz(1.3868706) q[3];
sx q[3];
rz(-1.9463332) q[3];
sx q[3];
rz(0.25115764) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-2.3544932) q[2];
sx q[2];
rz(-3.0566065) q[2];
sx q[2];
rz(-0.3698012) q[2];
rz(-1.9434816) q[3];
sx q[3];
rz(-1.5503649) q[3];
sx q[3];
rz(-0.91355356) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.13854606) q[0];
sx q[0];
rz(-1.9771165) q[0];
sx q[0];
rz(1.8856915) q[0];
rz(2.1239602) q[1];
sx q[1];
rz(-1.7435278) q[1];
sx q[1];
rz(1.7494038) q[1];
rz(3.0690003) q[2];
sx q[2];
rz(-1.6410927) q[2];
sx q[2];
rz(-2.1319364) q[2];
rz(0.16416488) q[3];
sx q[3];
rz(-1.4420684) q[3];
sx q[3];
rz(1.3185929) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
