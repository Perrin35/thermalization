OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-2.4797526) q[0];
sx q[0];
rz(-2.2979484) q[0];
sx q[0];
rz(2.9736829) q[0];
rz(1.1711988) q[1];
sx q[1];
rz(-2.8462703) q[1];
sx q[1];
rz(0.056161031) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0060318) q[0];
sx q[0];
rz(-2.6132085) q[0];
sx q[0];
rz(2.0023268) q[0];
rz(-0.21284717) q[2];
sx q[2];
rz(-2.2058862) q[2];
sx q[2];
rz(2.0095306) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.72370428) q[1];
sx q[1];
rz(-2.1065518) q[1];
sx q[1];
rz(2.8292311) q[1];
rz(-1.5264741) q[3];
sx q[3];
rz(-1.830415) q[3];
sx q[3];
rz(2.0093902) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(2.7636259) q[2];
sx q[2];
rz(-0.28181919) q[2];
sx q[2];
rz(-2.7089233) q[2];
rz(-1.1928605) q[3];
sx q[3];
rz(-1.2377219) q[3];
sx q[3];
rz(2.7584934) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.52779657) q[0];
sx q[0];
rz(-0.48848099) q[0];
sx q[0];
rz(1.8288076) q[0];
rz(0.20547543) q[1];
sx q[1];
rz(-2.165129) q[1];
sx q[1];
rz(-1.9899433) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5702471) q[0];
sx q[0];
rz(-0.76128188) q[0];
sx q[0];
rz(2.0061357) q[0];
rz(2.4685681) q[2];
sx q[2];
rz(-2.4403799) q[2];
sx q[2];
rz(-2.6075624) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(0.68576751) q[1];
sx q[1];
rz(-2.0058504) q[1];
sx q[1];
rz(1.1063834) q[1];
rz(-pi) q[2];
x q[2];
rz(0.63668164) q[3];
sx q[3];
rz(-2.5425306) q[3];
sx q[3];
rz(2.3707795) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-0.0097222086) q[2];
sx q[2];
rz(-1.4902318) q[2];
sx q[2];
rz(-2.9197664) q[2];
rz(-2.7644073) q[3];
sx q[3];
rz(-0.42755869) q[3];
sx q[3];
rz(0.77243531) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[3];
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
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8310228) q[0];
sx q[0];
rz(-3.0492058) q[0];
sx q[0];
rz(-0.036852766) q[0];
rz(0.82551461) q[1];
sx q[1];
rz(-1.3157536) q[1];
sx q[1];
rz(3.085014) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.53783137) q[0];
sx q[0];
rz(-2.1268401) q[0];
sx q[0];
rz(-2.6105196) q[0];
rz(-pi) q[1];
x q[1];
rz(2.668004) q[2];
sx q[2];
rz(-0.28652546) q[2];
sx q[2];
rz(2.5055656) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.0128855) q[1];
sx q[1];
rz(-1.371908) q[1];
sx q[1];
rz(-2.835564) q[1];
x q[2];
rz(-0.64485456) q[3];
sx q[3];
rz(-1.2554902) q[3];
sx q[3];
rz(-2.9063318) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-2.8686707) q[2];
sx q[2];
rz(-1.4977095) q[2];
sx q[2];
rz(0.92612129) q[2];
rz(-2.5849294) q[3];
sx q[3];
rz(-0.29354468) q[3];
sx q[3];
rz(2.0986957) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2264003) q[0];
sx q[0];
rz(-2.4214348) q[0];
sx q[0];
rz(-2.3994989) q[0];
rz(1.1391976) q[1];
sx q[1];
rz(-2.6622055) q[1];
sx q[1];
rz(2.6779968) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.050042) q[0];
sx q[0];
rz(-0.76489641) q[0];
sx q[0];
rz(2.0043623) q[0];
x q[1];
rz(-2.9651627) q[2];
sx q[2];
rz(-2.8039805) q[2];
sx q[2];
rz(0.38844973) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.5220118) q[1];
sx q[1];
rz(-1.7420235) q[1];
sx q[1];
rz(2.3492858) q[1];
x q[2];
rz(2.2673244) q[3];
sx q[3];
rz(-0.6797176) q[3];
sx q[3];
rz(3.0391039) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-1.7745557) q[2];
sx q[2];
rz(-0.15731263) q[2];
sx q[2];
rz(0.049499361) q[2];
rz(3.0130623) q[3];
sx q[3];
rz(-1.5481719) q[3];
sx q[3];
rz(-3.0226504) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0054935) q[0];
sx q[0];
rz(-0.43731421) q[0];
sx q[0];
rz(-2.8438925) q[0];
rz(0.4822576) q[1];
sx q[1];
rz(-2.3869956) q[1];
sx q[1];
rz(2.1972426) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0194861) q[0];
sx q[0];
rz(-1.5257611) q[0];
sx q[0];
rz(-1.7800063) q[0];
x q[1];
rz(-1.1733426) q[2];
sx q[2];
rz(-0.83565088) q[2];
sx q[2];
rz(3.1189001) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(0.22630616) q[1];
sx q[1];
rz(-1.5686791) q[1];
sx q[1];
rz(1.5096942) q[1];
rz(-pi) q[2];
rz(-1.827042) q[3];
sx q[3];
rz(-1.991193) q[3];
sx q[3];
rz(-2.8375569) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-1.258761) q[2];
sx q[2];
rz(-2.0488887) q[2];
sx q[2];
rz(3.0333701) q[2];
rz(0.0023068874) q[3];
sx q[3];
rz(-1.5294411) q[3];
sx q[3];
rz(-0.32430696) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
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
rz(-0.43679431) q[0];
sx q[0];
rz(-2.695485) q[0];
sx q[0];
rz(2.5571402) q[0];
rz(-2.2553518) q[1];
sx q[1];
rz(-0.61683547) q[1];
sx q[1];
rz(-0.054919682) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6829837) q[0];
sx q[0];
rz(-0.28124547) q[0];
sx q[0];
rz(-1.2693229) q[0];
x q[1];
rz(0.018718406) q[2];
sx q[2];
rz(-1.8939549) q[2];
sx q[2];
rz(1.2013555) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.37327787) q[1];
sx q[1];
rz(-1.0265961) q[1];
sx q[1];
rz(1.1276223) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.6782126) q[3];
sx q[3];
rz(-2.1043092) q[3];
sx q[3];
rz(-2.2802071) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-0.75446689) q[2];
sx q[2];
rz(-0.12067623) q[2];
sx q[2];
rz(-2.1248655) q[2];
rz(-0.54404849) q[3];
sx q[3];
rz(-0.35990158) q[3];
sx q[3];
rz(-1.8090766) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6761557) q[0];
sx q[0];
rz(-2.1570719) q[0];
sx q[0];
rz(-0.28453919) q[0];
rz(-0.94447213) q[1];
sx q[1];
rz(-1.9453134) q[1];
sx q[1];
rz(2.231266) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1658265) q[0];
sx q[0];
rz(-1.5162139) q[0];
sx q[0];
rz(-1.5058917) q[0];
rz(-pi) q[1];
rz(-0.34142999) q[2];
sx q[2];
rz(-1.176398) q[2];
sx q[2];
rz(0.1614801) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-1.4714204) q[1];
sx q[1];
rz(-0.52392611) q[1];
sx q[1];
rz(-0.36797221) q[1];
rz(-pi) q[2];
rz(0.99366412) q[3];
sx q[3];
rz(-0.93615195) q[3];
sx q[3];
rz(1.5821379) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.7515144) q[2];
sx q[2];
rz(-0.063515924) q[2];
sx q[2];
rz(0.92203036) q[2];
rz(0.56728029) q[3];
sx q[3];
rz(-1.4138979) q[3];
sx q[3];
rz(2.1218307) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.89408016) q[0];
sx q[0];
rz(-2.4649354) q[0];
sx q[0];
rz(-3.0122053) q[0];
rz(-0.63240504) q[1];
sx q[1];
rz(-2.1148732) q[1];
sx q[1];
rz(-0.30050373) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.520641) q[0];
sx q[0];
rz(-1.2004939) q[0];
sx q[0];
rz(0.83604367) q[0];
rz(0.95327611) q[2];
sx q[2];
rz(-0.91070181) q[2];
sx q[2];
rz(0.88027871) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-0.44725542) q[1];
sx q[1];
rz(-1.9112504) q[1];
sx q[1];
rz(1.7561595) q[1];
rz(1.237805) q[3];
sx q[3];
rz(-0.78562842) q[3];
sx q[3];
rz(1.9612519) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-0.58632103) q[2];
sx q[2];
rz(-2.1466612) q[2];
sx q[2];
rz(2.3596181) q[2];
rz(-0.55109763) q[3];
sx q[3];
rz(-1.3827773) q[3];
sx q[3];
rz(-0.15792318) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[2];
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
rz(-1.5732116) q[0];
sx q[0];
rz(-1.1567572) q[0];
sx q[0];
rz(-3.0138299) q[0];
rz(2.5993775) q[1];
sx q[1];
rz(-2.1844889) q[1];
sx q[1];
rz(-0.75884563) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1369143) q[0];
sx q[0];
rz(-1.5537964) q[0];
sx q[0];
rz(-1.1466115) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.3779638) q[2];
sx q[2];
rz(-2.170993) q[2];
sx q[2];
rz(-0.45229518) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-0.54770494) q[1];
sx q[1];
rz(-2.2624359) q[1];
sx q[1];
rz(1.3681075) q[1];
x q[2];
rz(3.127029) q[3];
sx q[3];
rz(-2.2363538) q[3];
sx q[3];
rz(1.2138106) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(1.1252497) q[2];
sx q[2];
rz(-1.3602076) q[2];
sx q[2];
rz(-0.49003595) q[2];
rz(1.4222493) q[3];
sx q[3];
rz(-1.2079206) q[3];
sx q[3];
rz(2.0786044) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
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
x q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.35995099) q[0];
sx q[0];
rz(-0.61976969) q[0];
sx q[0];
rz(0.075335659) q[0];
rz(0.8967337) q[1];
sx q[1];
rz(-1.1743841) q[1];
sx q[1];
rz(-0.60992253) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.9655351) q[0];
sx q[0];
rz(-0.81747222) q[0];
sx q[0];
rz(-0.39359351) q[0];
rz(-1.8462734) q[2];
sx q[2];
rz(-1.9343978) q[2];
sx q[2];
rz(1.7737349) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(0.016644195) q[1];
sx q[1];
rz(-2.2955107) q[1];
sx q[1];
rz(1.8925341) q[1];
rz(-pi) q[2];
rz(0.53819733) q[3];
sx q[3];
rz(-0.81848577) q[3];
sx q[3];
rz(1.324211) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-0.23218368) q[2];
sx q[2];
rz(-0.82513088) q[2];
sx q[2];
rz(-2.4278736) q[2];
rz(0.37832007) q[3];
sx q[3];
rz(-0.49351966) q[3];
sx q[3];
rz(2.2617214) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.80355766) q[0];
sx q[0];
rz(-1.150158) q[0];
sx q[0];
rz(-1.5858571) q[0];
rz(-2.4907885) q[1];
sx q[1];
rz(-1.6497859) q[1];
sx q[1];
rz(-0.12129687) q[1];
rz(-1.3748319) q[2];
sx q[2];
rz(-1.5407731) q[2];
sx q[2];
rz(-2.4827448) q[2];
rz(2.0402504) q[3];
sx q[3];
rz(-2.6228842) q[3];
sx q[3];
rz(1.7026671) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];