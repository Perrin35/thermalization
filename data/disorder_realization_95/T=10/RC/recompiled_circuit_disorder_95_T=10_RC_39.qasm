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
rz(5.2810623) q[0];
sx q[0];
rz(5.3856344) q[0];
rz(-3.3759723) q[1];
sx q[1];
rz(3.4174089) q[1];
sx q[1];
rz(13.630907) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9165186) q[0];
sx q[0];
rz(-1.6406035) q[0];
sx q[0];
rz(0.9401456) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.0374374) q[2];
sx q[2];
rz(-2.4541514) q[2];
sx q[2];
rz(-1.3165557) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-1.2805466) q[1];
sx q[1];
rz(-1.4006873) q[1];
sx q[1];
rz(-0.64330805) q[1];
rz(-pi) q[2];
rz(-1.7552056) q[3];
sx q[3];
rz(-0.82114906) q[3];
sx q[3];
rz(2.8781995) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-2.477318) q[2];
sx q[2];
rz(-0.6330108) q[2];
sx q[2];
rz(-0.73195362) q[2];
rz(-2.1814363) q[3];
sx q[3];
rz(-2.319016) q[3];
sx q[3];
rz(1.4320954) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.17523781) q[0];
sx q[0];
rz(-2.2580999) q[0];
sx q[0];
rz(0.91645855) q[0];
rz(0.48049277) q[1];
sx q[1];
rz(-2.5669211) q[1];
sx q[1];
rz(2.2629471) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6966349) q[0];
sx q[0];
rz(-0.96224552) q[0];
sx q[0];
rz(-2.6702325) q[0];
rz(2.5231045) q[2];
sx q[2];
rz(-1.3534708) q[2];
sx q[2];
rz(-2.6478812) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(0.52777427) q[1];
sx q[1];
rz(-1.6442181) q[1];
sx q[1];
rz(2.7223177) q[1];
rz(-pi) q[2];
rz(1.2611748) q[3];
sx q[3];
rz(-1.6358346) q[3];
sx q[3];
rz(-2.0464954) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-1.8615222) q[2];
sx q[2];
rz(-1.7333938) q[2];
sx q[2];
rz(-2.4943165) q[2];
rz(-0.17368008) q[3];
sx q[3];
rz(-2.0208385) q[3];
sx q[3];
rz(2.98996) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.52755255) q[0];
sx q[0];
rz(-0.63967597) q[0];
sx q[0];
rz(-0.24599427) q[0];
rz(-1.7315158) q[1];
sx q[1];
rz(-1.1743816) q[1];
sx q[1];
rz(-2.0203967) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4778053) q[0];
sx q[0];
rz(-2.2664547) q[0];
sx q[0];
rz(2.0545161) q[0];
x q[1];
rz(-0.063935117) q[2];
sx q[2];
rz(-1.9629515) q[2];
sx q[2];
rz(-2.8901697) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.1998636) q[1];
sx q[1];
rz(-0.45048303) q[1];
sx q[1];
rz(-2.4099039) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.6882012) q[3];
sx q[3];
rz(-2.061764) q[3];
sx q[3];
rz(-1.1743869) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-0.40144172) q[2];
sx q[2];
rz(-1.4462877) q[2];
sx q[2];
rz(-2.1901954) q[2];
rz(-2.4915063) q[3];
sx q[3];
rz(-1.254046) q[3];
sx q[3];
rz(0.29561177) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
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
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
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
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7096036) q[0];
sx q[0];
rz(-0.63392144) q[0];
sx q[0];
rz(3.0984127) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.94159796) q[2];
sx q[2];
rz(-2.5430352) q[2];
sx q[2];
rz(1.7184005) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(1.1981922) q[1];
sx q[1];
rz(-2.0640089) q[1];
sx q[1];
rz(2.0342159) q[1];
rz(-0.096480358) q[3];
sx q[3];
rz(-2.6498142) q[3];
sx q[3];
rz(1.552856) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-1.5300166) q[2];
sx q[2];
rz(-1.896603) q[2];
sx q[2];
rz(-2.8895565) q[2];
rz(-2.7633372) q[3];
sx q[3];
rz(-2.9791322) q[3];
sx q[3];
rz(-2.7799515) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.56548059) q[0];
sx q[0];
rz(-1.4030554) q[0];
sx q[0];
rz(0.53260032) q[0];
rz(-1.4415007) q[1];
sx q[1];
rz(-2.7756727) q[1];
sx q[1];
rz(1.9929569) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0557077) q[0];
sx q[0];
rz(-2.5143393) q[0];
sx q[0];
rz(1.6968615) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.2448036) q[2];
sx q[2];
rz(-1.2090346) q[2];
sx q[2];
rz(1.0587495) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(0.27302882) q[1];
sx q[1];
rz(-1.3820573) q[1];
sx q[1];
rz(0.36004685) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.1762189) q[3];
sx q[3];
rz(-0.88486457) q[3];
sx q[3];
rz(-1.8358313) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-2.1308412) q[2];
sx q[2];
rz(-0.66117078) q[2];
sx q[2];
rz(1.032069) q[2];
rz(0.71470913) q[3];
sx q[3];
rz(-1.2811477) q[3];
sx q[3];
rz(0.023795279) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.59153581) q[0];
sx q[0];
rz(-1.93601) q[0];
sx q[0];
rz(-0.39598879) q[0];
rz(1.4453325) q[1];
sx q[1];
rz(-1.4424125) q[1];
sx q[1];
rz(2.9352303) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9128742) q[0];
sx q[0];
rz(-1.657907) q[0];
sx q[0];
rz(2.3939783) q[0];
rz(-pi) q[1];
rz(2.2588737) q[2];
sx q[2];
rz(-2.1157017) q[2];
sx q[2];
rz(0.27872745) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.5620835) q[1];
sx q[1];
rz(-1.675203) q[1];
sx q[1];
rz(-3.0444006) q[1];
rz(-1.8984853) q[3];
sx q[3];
rz(-1.7565691) q[3];
sx q[3];
rz(0.9542619) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(0.53987327) q[2];
sx q[2];
rz(-1.0605992) q[2];
sx q[2];
rz(-2.8708141) q[2];
rz(0.74832908) q[3];
sx q[3];
rz(-1.3724519) q[3];
sx q[3];
rz(1.07871) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8966184) q[0];
sx q[0];
rz(-2.0386319) q[0];
sx q[0];
rz(0.57624972) q[0];
rz(0.17310625) q[1];
sx q[1];
rz(-0.74179596) q[1];
sx q[1];
rz(-1.2111838) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1887218) q[0];
sx q[0];
rz(-2.0824008) q[0];
sx q[0];
rz(-2.1923724) q[0];
rz(-1.7089825) q[2];
sx q[2];
rz(-2.0565363) q[2];
sx q[2];
rz(-1.9463584) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-2.3019575) q[1];
sx q[1];
rz(-2.9710899) q[1];
sx q[1];
rz(1.9819928) q[1];
rz(-pi) q[2];
rz(0.96774775) q[3];
sx q[3];
rz(-1.3021886) q[3];
sx q[3];
rz(-0.3097765) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(1.8048191) q[2];
sx q[2];
rz(-2.4833198) q[2];
sx q[2];
rz(-3.1170735) q[2];
rz(0.71497861) q[3];
sx q[3];
rz(-1.4638126) q[3];
sx q[3];
rz(1.5589176) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
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
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.0054758469) q[0];
sx q[0];
rz(-1.550721) q[0];
sx q[0];
rz(-2.1210282) q[0];
rz(0.15469805) q[1];
sx q[1];
rz(-1.5203412) q[1];
sx q[1];
rz(-1.9205836) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.3295791) q[0];
sx q[0];
rz(-2.3276969) q[0];
sx q[0];
rz(-1.194186) q[0];
rz(-0.17692716) q[2];
sx q[2];
rz(-1.7311586) q[2];
sx q[2];
rz(-0.71080506) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(0.11215969) q[1];
sx q[1];
rz(-2.9661848) q[1];
sx q[1];
rz(0.68316858) q[1];
rz(2.3200795) q[3];
sx q[3];
rz(-1.320208) q[3];
sx q[3];
rz(0.65419765) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-2.8528379) q[2];
sx q[2];
rz(-1.6483665) q[2];
sx q[2];
rz(-0.70518804) q[2];
rz(0.60338902) q[3];
sx q[3];
rz(-2.2796977) q[3];
sx q[3];
rz(-1.5195297) q[3];
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
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0309546) q[0];
sx q[0];
rz(-1.5752666) q[0];
sx q[0];
rz(-3.0045793) q[0];
rz(-0.6048454) q[1];
sx q[1];
rz(-0.73692656) q[1];
sx q[1];
rz(-3.0158214) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9916358) q[0];
sx q[0];
rz(-1.7466674) q[0];
sx q[0];
rz(-1.3791023) q[0];
x q[1];
rz(-1.7681098) q[2];
sx q[2];
rz(-1.8981877) q[2];
sx q[2];
rz(-1.2158074) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(2.3667664) q[1];
sx q[1];
rz(-0.80087304) q[1];
sx q[1];
rz(-1.1714539) q[1];
rz(-pi) q[2];
x q[2];
rz(0.59154193) q[3];
sx q[3];
rz(-0.35174832) q[3];
sx q[3];
rz(1.8464309) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(3.003309) q[2];
sx q[2];
rz(-0.97390276) q[2];
sx q[2];
rz(-2.4323145) q[2];
rz(0.62018958) q[3];
sx q[3];
rz(-0.87738335) q[3];
sx q[3];
rz(0.15670776) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9897292) q[0];
sx q[0];
rz(-0.79280889) q[0];
sx q[0];
rz(1.9412769) q[0];
rz(-2.8740846) q[1];
sx q[1];
rz(-2.2876883) q[1];
sx q[1];
rz(2.0013924) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8126497) q[0];
sx q[0];
rz(-0.2812627) q[0];
sx q[0];
rz(1.2055231) q[0];
x q[1];
rz(0.26009772) q[2];
sx q[2];
rz(-2.0195228) q[2];
sx q[2];
rz(1.5425494) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-0.89123594) q[1];
sx q[1];
rz(-1.48238) q[1];
sx q[1];
rz(0.18914117) q[1];
rz(-pi) q[2];
x q[2];
rz(2.8285698) q[3];
sx q[3];
rz(-0.81837294) q[3];
sx q[3];
rz(-0.5648053) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-2.7662979) q[2];
sx q[2];
rz(-1.6833498) q[2];
sx q[2];
rz(0.89861384) q[2];
rz(0.0071772655) q[3];
sx q[3];
rz(-2.3982748) q[3];
sx q[3];
rz(1.7858508) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.89467775) q[0];
sx q[0];
rz(-1.7264195) q[0];
sx q[0];
rz(1.3608426) q[0];
rz(0.5207516) q[1];
sx q[1];
rz(-1.3856577) q[1];
sx q[1];
rz(1.7210977) q[1];
rz(-1.0976787) q[2];
sx q[2];
rz(-1.1259148) q[2];
sx q[2];
rz(-1.7150707) q[2];
rz(-2.5614212) q[3];
sx q[3];
rz(-0.67065722) q[3];
sx q[3];
rz(-1.8967659) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
