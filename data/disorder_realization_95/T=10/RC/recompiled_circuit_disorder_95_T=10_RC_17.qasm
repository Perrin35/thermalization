OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-0.99958217) q[0];
sx q[0];
rz(-1.002123) q[0];
sx q[0];
rz(2.2440417) q[0];
rz(-0.23437962) q[1];
sx q[1];
rz(-0.27581629) q[1];
sx q[1];
rz(2.0770567) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.25047725) q[0];
sx q[0];
rz(-2.5076137) q[0];
sx q[0];
rz(-1.4527713) q[0];
rz(-pi) q[1];
rz(-1.1041553) q[2];
sx q[2];
rz(-0.68744126) q[2];
sx q[2];
rz(1.825037) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-1.2805466) q[1];
sx q[1];
rz(-1.4006873) q[1];
sx q[1];
rz(2.4982846) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.7552056) q[3];
sx q[3];
rz(-2.3204436) q[3];
sx q[3];
rz(0.26339312) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(0.66427461) q[2];
sx q[2];
rz(-2.5085818) q[2];
sx q[2];
rz(0.73195362) q[2];
rz(-0.96015635) q[3];
sx q[3];
rz(-0.82257661) q[3];
sx q[3];
rz(-1.7094973) q[3];
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
sx q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.17523781) q[0];
sx q[0];
rz(-0.8834928) q[0];
sx q[0];
rz(0.91645855) q[0];
rz(-0.48049277) q[1];
sx q[1];
rz(-2.5669211) q[1];
sx q[1];
rz(-2.2629471) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.71683305) q[0];
sx q[0];
rz(-0.75100198) q[0];
sx q[0];
rz(-0.99320937) q[0];
x q[1];
rz(-1.3061366) q[2];
sx q[2];
rz(-0.96894962) q[2];
sx q[2];
rz(2.2167609) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-0.87989992) q[1];
sx q[1];
rz(-0.42527929) q[1];
sx q[1];
rz(0.17875032) q[1];
rz(-1.3602123) q[3];
sx q[3];
rz(-2.8254291) q[3];
sx q[3];
rz(-2.8663243) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-1.2800704) q[2];
sx q[2];
rz(-1.7333938) q[2];
sx q[2];
rz(-0.64727616) q[2];
rz(0.17368008) q[3];
sx q[3];
rz(-1.1207542) q[3];
sx q[3];
rz(-0.15163264) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6140401) q[0];
sx q[0];
rz(-0.63967597) q[0];
sx q[0];
rz(0.24599427) q[0];
rz(-1.4100769) q[1];
sx q[1];
rz(-1.9672111) q[1];
sx q[1];
rz(1.1211959) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4778053) q[0];
sx q[0];
rz(-2.2664547) q[0];
sx q[0];
rz(1.0870766) q[0];
rz(-pi) q[1];
rz(-1.4175225) q[2];
sx q[2];
rz(-0.3970662) q[2];
sx q[2];
rz(3.0561471) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-1.9417291) q[1];
sx q[1];
rz(-2.6911096) q[1];
sx q[1];
rz(-0.73168879) q[1];
x q[2];
rz(-2.9259053) q[3];
sx q[3];
rz(-0.50369278) q[3];
sx q[3];
rz(-1.7220725) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-0.40144172) q[2];
sx q[2];
rz(-1.4462877) q[2];
sx q[2];
rz(0.95139727) q[2];
rz(2.4915063) q[3];
sx q[3];
rz(-1.8875467) q[3];
sx q[3];
rz(0.29561177) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.51405108) q[0];
sx q[0];
rz(-2.5294332) q[0];
sx q[0];
rz(2.3535368) q[0];
rz(1.0568985) q[1];
sx q[1];
rz(-1.1833271) q[1];
sx q[1];
rz(-3.025211) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.43198904) q[0];
sx q[0];
rz(-2.5076712) q[0];
sx q[0];
rz(-3.0984127) q[0];
x q[1];
rz(2.0747244) q[2];
sx q[2];
rz(-1.908784) q[2];
sx q[2];
rz(-2.4525814) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-1.9434005) q[1];
sx q[1];
rz(-2.0640089) q[1];
sx q[1];
rz(2.0342159) q[1];
rz(-pi) q[2];
x q[2];
rz(-3.0451123) q[3];
sx q[3];
rz(-2.6498142) q[3];
sx q[3];
rz(-1.552856) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.6115761) q[2];
sx q[2];
rz(-1.2449896) q[2];
sx q[2];
rz(2.8895565) q[2];
rz(-0.37825545) q[3];
sx q[3];
rz(-0.16246048) q[3];
sx q[3];
rz(0.3616412) q[3];
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
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.56548059) q[0];
sx q[0];
rz(-1.7385372) q[0];
sx q[0];
rz(2.6089923) q[0];
rz(1.700092) q[1];
sx q[1];
rz(-0.36591995) q[1];
sx q[1];
rz(-1.9929569) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.085885) q[0];
sx q[0];
rz(-2.5143393) q[0];
sx q[0];
rz(1.6968615) q[0];
rz(2.2448036) q[2];
sx q[2];
rz(-1.9325581) q[2];
sx q[2];
rz(1.0587495) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-2.8685638) q[1];
sx q[1];
rz(-1.7595353) q[1];
sx q[1];
rz(-0.36004685) q[1];
rz(-pi) q[2];
rz(2.4162021) q[3];
sx q[3];
rz(-1.2687506) q[3];
sx q[3];
rz(-0.52291742) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-1.0107515) q[2];
sx q[2];
rz(-2.4804219) q[2];
sx q[2];
rz(-2.1095236) q[2];
rz(-2.4268835) q[3];
sx q[3];
rz(-1.8604449) q[3];
sx q[3];
rz(3.1177974) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
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
rz(-pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.59153581) q[0];
sx q[0];
rz(-1.2055826) q[0];
sx q[0];
rz(-2.7456039) q[0];
rz(-1.4453325) q[1];
sx q[1];
rz(-1.6991801) q[1];
sx q[1];
rz(2.9352303) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8800031) q[0];
sx q[0];
rz(-0.82669175) q[0];
sx q[0];
rz(-1.4522626) q[0];
rz(-pi) q[1];
x q[1];
rz(0.80870734) q[2];
sx q[2];
rz(-2.2924097) q[2];
sx q[2];
rz(-1.2869814) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-2.5620835) q[1];
sx q[1];
rz(-1.675203) q[1];
sx q[1];
rz(-3.0444006) q[1];
rz(1.0422802) q[3];
sx q[3];
rz(-2.7665666) q[3];
sx q[3];
rz(3.022775) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(2.6017194) q[2];
sx q[2];
rz(-1.0605992) q[2];
sx q[2];
rz(2.8708141) q[2];
rz(0.74832908) q[3];
sx q[3];
rz(-1.3724519) q[3];
sx q[3];
rz(-2.0628827) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8966184) q[0];
sx q[0];
rz(-2.0386319) q[0];
sx q[0];
rz(-2.5653429) q[0];
rz(2.9684864) q[1];
sx q[1];
rz(-2.3997967) q[1];
sx q[1];
rz(1.9304088) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9528708) q[0];
sx q[0];
rz(-2.0824008) q[0];
sx q[0];
rz(-2.1923724) q[0];
rz(-pi) q[1];
rz(-2.8863635) q[2];
sx q[2];
rz(-0.50349871) q[2];
sx q[2];
rz(2.2358759) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-2.8163029) q[1];
sx q[1];
rz(-1.6386697) q[1];
sx q[1];
rz(-1.7273278) q[1];
x q[2];
rz(-2.8190523) q[3];
sx q[3];
rz(-2.1493559) q[3];
sx q[3];
rz(-1.4417779) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(1.3367735) q[2];
sx q[2];
rz(-2.4833198) q[2];
sx q[2];
rz(3.1170735) q[2];
rz(-0.71497861) q[3];
sx q[3];
rz(-1.67778) q[3];
sx q[3];
rz(1.5589176) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.1361168) q[0];
sx q[0];
rz(-1.550721) q[0];
sx q[0];
rz(1.0205644) q[0];
rz(-0.15469805) q[1];
sx q[1];
rz(-1.6212515) q[1];
sx q[1];
rz(-1.9205836) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.97604254) q[0];
sx q[0];
rz(-1.8414458) q[0];
sx q[0];
rz(2.3483777) q[0];
rz(-pi) q[1];
rz(0.74322015) q[2];
sx q[2];
rz(-0.23822242) q[2];
sx q[2];
rz(1.5889578) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-3.029433) q[1];
sx q[1];
rz(-2.9661848) q[1];
sx q[1];
rz(-2.4584241) q[1];
rz(0.82151316) q[3];
sx q[3];
rz(-1.320208) q[3];
sx q[3];
rz(2.487395) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(2.8528379) q[2];
sx q[2];
rz(-1.4932262) q[2];
sx q[2];
rz(2.4364046) q[2];
rz(-0.60338902) q[3];
sx q[3];
rz(-0.861895) q[3];
sx q[3];
rz(1.6220629) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(pi/2) q[3];
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
rz(-2.0309546) q[0];
sx q[0];
rz(-1.5663261) q[0];
sx q[0];
rz(0.13701339) q[0];
rz(0.6048454) q[1];
sx q[1];
rz(-2.4046661) q[1];
sx q[1];
rz(0.12577122) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.754697) q[0];
sx q[0];
rz(-1.3820952) q[0];
sx q[0];
rz(0.17908355) q[0];
rz(-0.52351064) q[2];
sx q[2];
rz(-0.38041174) q[2];
sx q[2];
rz(-1.7720122) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-2.3667664) q[1];
sx q[1];
rz(-2.3407196) q[1];
sx q[1];
rz(1.9701387) q[1];
rz(-pi) q[2];
x q[2];
rz(2.5500507) q[3];
sx q[3];
rz(-0.35174832) q[3];
sx q[3];
rz(-1.8464309) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(0.13828364) q[2];
sx q[2];
rz(-0.97390276) q[2];
sx q[2];
rz(2.4323145) q[2];
rz(-0.62018958) q[3];
sx q[3];
rz(-0.87738335) q[3];
sx q[3];
rz(-0.15670776) q[3];
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
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1518635) q[0];
sx q[0];
rz(-2.3487838) q[0];
sx q[0];
rz(-1.2003157) q[0];
rz(-0.26750803) q[1];
sx q[1];
rz(-2.2876883) q[1];
sx q[1];
rz(1.1402003) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.32894293) q[0];
sx q[0];
rz(-0.2812627) q[0];
sx q[0];
rz(-1.2055231) q[0];
x q[1];
rz(0.26009772) q[2];
sx q[2];
rz(-1.1220699) q[2];
sx q[2];
rz(1.5990433) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-0.89123594) q[1];
sx q[1];
rz(-1.48238) q[1];
sx q[1];
rz(-2.9524515) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.2530008) q[3];
sx q[3];
rz(-0.80298775) q[3];
sx q[3];
rz(-0.12249891) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.37529477) q[2];
sx q[2];
rz(-1.4582429) q[2];
sx q[2];
rz(-2.2429788) q[2];
rz(-3.1344154) q[3];
sx q[3];
rz(-0.74331784) q[3];
sx q[3];
rz(1.3557419) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2469149) q[0];
sx q[0];
rz(-1.7264195) q[0];
sx q[0];
rz(1.3608426) q[0];
rz(-2.6208411) q[1];
sx q[1];
rz(-1.3856577) q[1];
sx q[1];
rz(1.7210977) q[1];
rz(2.3788135) q[2];
sx q[2];
rz(-0.63763466) q[2];
sx q[2];
rz(0.55479738) q[2];
rz(-0.58017147) q[3];
sx q[3];
rz(-2.4709354) q[3];
sx q[3];
rz(1.2448268) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
