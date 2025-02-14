OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(0.98915339) q[0];
sx q[0];
rz(1.5751155) q[0];
sx q[0];
rz(10.549904) q[0];
rz(-2.1195124) q[1];
sx q[1];
rz(-2.4740969) q[1];
sx q[1];
rz(-0.8134841) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.33685499) q[0];
sx q[0];
rz(-1.2053524) q[0];
sx q[0];
rz(-2.9904537) q[0];
x q[1];
rz(1.9314693) q[2];
sx q[2];
rz(-1.9934335) q[2];
sx q[2];
rz(0.3061184) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(0.58373815) q[1];
sx q[1];
rz(-2.0101133) q[1];
sx q[1];
rz(1.735461) q[1];
rz(-1.068068) q[3];
sx q[3];
rz(-2.9554071) q[3];
sx q[3];
rz(-1.26621) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(0.4890613) q[2];
sx q[2];
rz(-1.6901313) q[2];
sx q[2];
rz(2.2129464) q[2];
rz(-1.5993902) q[3];
sx q[3];
rz(-1.8079115) q[3];
sx q[3];
rz(1.1716051) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9376675) q[0];
sx q[0];
rz(-1.7610022) q[0];
sx q[0];
rz(0.12705886) q[0];
rz(2.1584885) q[1];
sx q[1];
rz(-1.3652722) q[1];
sx q[1];
rz(2.3703221) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.49202) q[0];
sx q[0];
rz(-0.82081534) q[0];
sx q[0];
rz(1.0787021) q[0];
rz(-1.2331687) q[2];
sx q[2];
rz(-1.513483) q[2];
sx q[2];
rz(-1.415451) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.018056695) q[1];
sx q[1];
rz(-0.50790826) q[1];
sx q[1];
rz(-2.3695709) q[1];
rz(-pi) q[2];
x q[2];
rz(2.3985474) q[3];
sx q[3];
rz(-0.72297572) q[3];
sx q[3];
rz(-1.522775) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-1.9942921) q[2];
sx q[2];
rz(-1.9376126) q[2];
sx q[2];
rz(-3.1331983) q[2];
rz(0.66347915) q[3];
sx q[3];
rz(-1.9050262) q[3];
sx q[3];
rz(0.26507637) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0405149) q[0];
sx q[0];
rz(-0.82656693) q[0];
sx q[0];
rz(-0.44152942) q[0];
rz(2.1532374) q[1];
sx q[1];
rz(-2.007273) q[1];
sx q[1];
rz(-3.006014) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.044047318) q[0];
sx q[0];
rz(-1.5545465) q[0];
sx q[0];
rz(1.5784997) q[0];
x q[1];
rz(-0.43164469) q[2];
sx q[2];
rz(-2.3767391) q[2];
sx q[2];
rz(-0.60827309) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-1.5244658) q[1];
sx q[1];
rz(-1.9248665) q[1];
sx q[1];
rz(-1.0431784) q[1];
x q[2];
rz(2.9445678) q[3];
sx q[3];
rz(-2.2954515) q[3];
sx q[3];
rz(2.0565095) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(2.2077937) q[2];
sx q[2];
rz(-1.4321045) q[2];
sx q[2];
rz(3.1033893) q[2];
rz(0.52538747) q[3];
sx q[3];
rz(-2.5140258) q[3];
sx q[3];
rz(2.4562522) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4811089) q[0];
sx q[0];
rz(-2.3479192) q[0];
sx q[0];
rz(-2.5307181) q[0];
rz(1.8065709) q[1];
sx q[1];
rz(-1.7673312) q[1];
sx q[1];
rz(-2.9023721) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6253875) q[0];
sx q[0];
rz(-2.5175178) q[0];
sx q[0];
rz(2.027987) q[0];
rz(-pi) q[1];
rz(-1.9630505) q[2];
sx q[2];
rz(-1.129727) q[2];
sx q[2];
rz(2.6792996) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-1.4799616) q[1];
sx q[1];
rz(-2.0212272) q[1];
sx q[1];
rz(-1.9305139) q[1];
rz(-1.7710822) q[3];
sx q[3];
rz(-1.1173751) q[3];
sx q[3];
rz(3.0016921) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(1.4227582) q[2];
sx q[2];
rz(-1.4449291) q[2];
sx q[2];
rz(-0.0017496721) q[2];
rz(-0.29378978) q[3];
sx q[3];
rz(-1.1894476) q[3];
sx q[3];
rz(-0.26688117) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9002429) q[0];
sx q[0];
rz(-1.7647864) q[0];
sx q[0];
rz(-2.2985261) q[0];
rz(-0.33755606) q[1];
sx q[1];
rz(-1.2755716) q[1];
sx q[1];
rz(-1.6036124) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0564954) q[0];
sx q[0];
rz(-0.81696327) q[0];
sx q[0];
rz(-2.8780065) q[0];
x q[1];
rz(-0.17436738) q[2];
sx q[2];
rz(-1.5182759) q[2];
sx q[2];
rz(2.8157521) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-1.5462436) q[1];
sx q[1];
rz(-1.3145295) q[1];
sx q[1];
rz(1.7728189) q[1];
rz(-pi) q[2];
x q[2];
rz(2.7790623) q[3];
sx q[3];
rz(-0.22386079) q[3];
sx q[3];
rz(0.26765841) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(2.7097077) q[2];
sx q[2];
rz(-1.0872492) q[2];
sx q[2];
rz(1.0464) q[2];
rz(2.8228068) q[3];
sx q[3];
rz(-2.4304978) q[3];
sx q[3];
rz(-0.54767245) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.043561291) q[0];
sx q[0];
rz(-1.8836319) q[0];
sx q[0];
rz(0.64796722) q[0];
rz(1.3994392) q[1];
sx q[1];
rz(-1.49767) q[1];
sx q[1];
rz(1.3290149) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.15547046) q[0];
sx q[0];
rz(-2.8948445) q[0];
sx q[0];
rz(-2.0339436) q[0];
rz(-1.7850661) q[2];
sx q[2];
rz(-1.2870711) q[2];
sx q[2];
rz(1.2905215) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-1.4635515) q[1];
sx q[1];
rz(-1.2604509) q[1];
sx q[1];
rz(-0.84872679) q[1];
rz(-pi) q[2];
x q[2];
rz(0.59698128) q[3];
sx q[3];
rz(-0.59904811) q[3];
sx q[3];
rz(2.0839276) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-0.38137388) q[2];
sx q[2];
rz(-1.9332644) q[2];
sx q[2];
rz(-1.5927429) q[2];
rz(-0.069843944) q[3];
sx q[3];
rz(-1.8826238) q[3];
sx q[3];
rz(0.86863345) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1182564) q[0];
sx q[0];
rz(-1.2160439) q[0];
sx q[0];
rz(2.1997531) q[0];
rz(-2.9815004) q[1];
sx q[1];
rz(-1.4804877) q[1];
sx q[1];
rz(-0.19217415) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.20142444) q[0];
sx q[0];
rz(-1.3469101) q[0];
sx q[0];
rz(-3.1046449) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.0490369) q[2];
sx q[2];
rz(-1.6055487) q[2];
sx q[2];
rz(3.1086189) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-0.46017473) q[1];
sx q[1];
rz(-0.81528403) q[1];
sx q[1];
rz(-2.6667206) q[1];
x q[2];
rz(-0.022054733) q[3];
sx q[3];
rz(-1.0733177) q[3];
sx q[3];
rz(2.9201404) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-2.3840702) q[2];
sx q[2];
rz(-1.9033868) q[2];
sx q[2];
rz(-2.7080217) q[2];
rz(1.5571669) q[3];
sx q[3];
rz(-1.4945364) q[3];
sx q[3];
rz(-2.7092194) q[3];
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
rz(-pi) q[0];
sx q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.50345355) q[0];
sx q[0];
rz(-1.3766377) q[0];
sx q[0];
rz(1.3442511) q[0];
rz(-1.9380219) q[1];
sx q[1];
rz(-1.9529587) q[1];
sx q[1];
rz(1.7787836) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8289611) q[0];
sx q[0];
rz(-1.5927218) q[0];
sx q[0];
rz(1.5410822) q[0];
rz(-1.7596471) q[2];
sx q[2];
rz(-0.7536234) q[2];
sx q[2];
rz(-0.058503956) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-2.2143933) q[1];
sx q[1];
rz(-1.7991369) q[1];
sx q[1];
rz(-1.6990927) q[1];
rz(-2.8239488) q[3];
sx q[3];
rz(-0.90932019) q[3];
sx q[3];
rz(2.3931488) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(0.56889304) q[2];
sx q[2];
rz(-2.3685444) q[2];
sx q[2];
rz(0.9838689) q[2];
rz(0.24108663) q[3];
sx q[3];
rz(-0.9328931) q[3];
sx q[3];
rz(-0.69055313) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6702061) q[0];
sx q[0];
rz(-0.79600483) q[0];
sx q[0];
rz(2.8669226) q[0];
rz(-1.341691) q[1];
sx q[1];
rz(-2.5119753) q[1];
sx q[1];
rz(1.0460269) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2316596) q[0];
sx q[0];
rz(-2.446736) q[0];
sx q[0];
rz(-0.52082638) q[0];
rz(-pi) q[1];
rz(-1.9555904) q[2];
sx q[2];
rz(-0.21285393) q[2];
sx q[2];
rz(1.720088) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-2.6089692) q[1];
sx q[1];
rz(-1.5057505) q[1];
sx q[1];
rz(2.5576127) q[1];
x q[2];
rz(1.0677797) q[3];
sx q[3];
rz(-1.5340337) q[3];
sx q[3];
rz(3.0779148) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(0.86299738) q[2];
sx q[2];
rz(-1.3906761) q[2];
sx q[2];
rz(0.39946237) q[2];
rz(0.56600371) q[3];
sx q[3];
rz(-2.1750906) q[3];
sx q[3];
rz(-0.81542265) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
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
rz(0.28867662) q[0];
sx q[0];
rz(-1.7221907) q[0];
sx q[0];
rz(-2.5592819) q[0];
rz(2.9583926) q[1];
sx q[1];
rz(-2.2661426) q[1];
sx q[1];
rz(0.80642548) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9574653) q[0];
sx q[0];
rz(-3.0054207) q[0];
sx q[0];
rz(0.86958142) q[0];
x q[1];
rz(-1.8381133) q[2];
sx q[2];
rz(-1.8952994) q[2];
sx q[2];
rz(2.7410206) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-0.90102531) q[1];
sx q[1];
rz(-1.7738924) q[1];
sx q[1];
rz(2.3939449) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.29911689) q[3];
sx q[3];
rz(-2.7795305) q[3];
sx q[3];
rz(2.5821834) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(2.9091984) q[2];
sx q[2];
rz(-2.4269203) q[2];
sx q[2];
rz(-0.8052899) q[2];
rz(-0.14690873) q[3];
sx q[3];
rz(-2.3555136) q[3];
sx q[3];
rz(0.73827353) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.88937) q[0];
sx q[0];
rz(-1.7563553) q[0];
sx q[0];
rz(1.8808543) q[0];
rz(-1.5421142) q[1];
sx q[1];
rz(-1.5972932) q[1];
sx q[1];
rz(1.6398026) q[1];
rz(-1.7003822) q[2];
sx q[2];
rz(-2.414188) q[2];
sx q[2];
rz(1.3368171) q[2];
rz(-0.44381683) q[3];
sx q[3];
rz(-2.4469821) q[3];
sx q[3];
rz(-1.2367005) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
