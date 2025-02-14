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
rz(0.17231365) q[0];
sx q[0];
rz(3.3495164) q[0];
sx q[0];
rz(10.715545) q[0];
rz(-2.9895904) q[1];
sx q[1];
rz(-0.64270371) q[1];
sx q[1];
rz(2.7132577) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1550395) q[0];
sx q[0];
rz(-1.6800092) q[0];
sx q[0];
rz(-3.0364353) q[0];
rz(-pi) q[1];
rz(-0.67240388) q[2];
sx q[2];
rz(-1.8032296) q[2];
sx q[2];
rz(1.642579) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(0.92378919) q[1];
sx q[1];
rz(-2.3084967) q[1];
sx q[1];
rz(-1.5125201) q[1];
rz(-2.9228404) q[3];
sx q[3];
rz(-1.7667337) q[3];
sx q[3];
rz(2.4234555) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(0.20994818) q[2];
sx q[2];
rz(-0.95982176) q[2];
sx q[2];
rz(1.2098562) q[2];
rz(3.0620388) q[3];
sx q[3];
rz(-2.7456561) q[3];
sx q[3];
rz(-2.615926) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[2];
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
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2570268) q[0];
sx q[0];
rz(-2.631812) q[0];
sx q[0];
rz(0.38401815) q[0];
rz(2.2513023) q[1];
sx q[1];
rz(-1.9046116) q[1];
sx q[1];
rz(0.30337897) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3502894) q[0];
sx q[0];
rz(-1.7614189) q[0];
sx q[0];
rz(2.1146569) q[0];
x q[1];
rz(2.0069507) q[2];
sx q[2];
rz(-0.56098962) q[2];
sx q[2];
rz(0.84535384) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.56306707) q[1];
sx q[1];
rz(-1.6614698) q[1];
sx q[1];
rz(3.0791326) q[1];
rz(-pi) q[2];
x q[2];
rz(2.4879913) q[3];
sx q[3];
rz(-1.2051688) q[3];
sx q[3];
rz(2.0662243) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(2.7552135) q[2];
sx q[2];
rz(-2.0342125) q[2];
sx q[2];
rz(2.4066822) q[2];
rz(-2.2928061) q[3];
sx q[3];
rz(-2.3990192) q[3];
sx q[3];
rz(-0.36062226) q[3];
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
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.81247771) q[0];
sx q[0];
rz(-2.766093) q[0];
sx q[0];
rz(-0.078711674) q[0];
rz(2.3178237) q[1];
sx q[1];
rz(-2.0562101) q[1];
sx q[1];
rz(-1.8471921) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1549653) q[0];
sx q[0];
rz(-1.6324348) q[0];
sx q[0];
rz(3.0810647) q[0];
x q[1];
rz(2.6316775) q[2];
sx q[2];
rz(-1.0140061) q[2];
sx q[2];
rz(-2.315588) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.8405806) q[1];
sx q[1];
rz(-0.30547474) q[1];
sx q[1];
rz(2.7423647) q[1];
x q[2];
rz(0.23842509) q[3];
sx q[3];
rz(-0.88078558) q[3];
sx q[3];
rz(-0.7499786) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-0.89982975) q[2];
sx q[2];
rz(-0.37909847) q[2];
sx q[2];
rz(-2.3251593) q[2];
rz(1.624931) q[3];
sx q[3];
rz(-2.1818325) q[3];
sx q[3];
rz(2.4412156) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.65275943) q[0];
sx q[0];
rz(-2.3214898) q[0];
sx q[0];
rz(-0.56831992) q[0];
rz(-1.142451) q[1];
sx q[1];
rz(-1.7894952) q[1];
sx q[1];
rz(-1.4521339) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.88599151) q[0];
sx q[0];
rz(-1.3491677) q[0];
sx q[0];
rz(2.5021573) q[0];
rz(-pi) q[1];
x q[1];
rz(1.4771039) q[2];
sx q[2];
rz(-1.3917275) q[2];
sx q[2];
rz(2.8833517) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(2.8206222) q[1];
sx q[1];
rz(-1.2577211) q[1];
sx q[1];
rz(-2.2105107) q[1];
rz(-pi) q[2];
rz(1.9098572) q[3];
sx q[3];
rz(-1.7393469) q[3];
sx q[3];
rz(-1.3685365) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-2.6292754) q[2];
sx q[2];
rz(-2.6111111) q[2];
sx q[2];
rz(0.076889195) q[2];
rz(2.7422089) q[3];
sx q[3];
rz(-0.9136343) q[3];
sx q[3];
rz(-0.54401773) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9870616) q[0];
sx q[0];
rz(-2.7847544) q[0];
sx q[0];
rz(0.29990184) q[0];
rz(-1.9918282) q[1];
sx q[1];
rz(-2.1297784) q[1];
sx q[1];
rz(-2.6531175) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.66118427) q[0];
sx q[0];
rz(-1.5380435) q[0];
sx q[0];
rz(1.7561654) q[0];
rz(-pi) q[1];
x q[1];
rz(2.0769507) q[2];
sx q[2];
rz(-2.1035668) q[2];
sx q[2];
rz(1.8440994) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(0.2327959) q[1];
sx q[1];
rz(-2.1543062) q[1];
sx q[1];
rz(1.5934492) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.4998922) q[3];
sx q[3];
rz(-1.3535168) q[3];
sx q[3];
rz(3.1171796) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-0.38146314) q[2];
sx q[2];
rz(-3.0035786) q[2];
sx q[2];
rz(-1.6628954) q[2];
rz(2.0315157) q[3];
sx q[3];
rz(-2.1434982) q[3];
sx q[3];
rz(2.4983675) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4019302) q[0];
sx q[0];
rz(-1.9597541) q[0];
sx q[0];
rz(-2.8231743) q[0];
rz(0.034612522) q[1];
sx q[1];
rz(-0.56762677) q[1];
sx q[1];
rz(-0.45077032) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0552669) q[0];
sx q[0];
rz(-2.0811715) q[0];
sx q[0];
rz(0.016252131) q[0];
x q[1];
rz(-2.320596) q[2];
sx q[2];
rz(-0.85047532) q[2];
sx q[2];
rz(-0.79325097) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(0.26397005) q[1];
sx q[1];
rz(-1.1671986) q[1];
sx q[1];
rz(2.2303459) q[1];
x q[2];
rz(1.1799906) q[3];
sx q[3];
rz(-2.3473661) q[3];
sx q[3];
rz(0.022217928) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-3.0938809) q[2];
sx q[2];
rz(-0.1897976) q[2];
sx q[2];
rz(-0.06165687) q[2];
rz(-3.0430074) q[3];
sx q[3];
rz(-0.74461377) q[3];
sx q[3];
rz(-2.4390167) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3848569) q[0];
sx q[0];
rz(-2.0282133) q[0];
sx q[0];
rz(3.1016438) q[0];
rz(0.7705676) q[1];
sx q[1];
rz(-2.4295085) q[1];
sx q[1];
rz(0.22824731) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4994482) q[0];
sx q[0];
rz(-3.0239883) q[0];
sx q[0];
rz(-2.2884877) q[0];
x q[1];
rz(1.2646152) q[2];
sx q[2];
rz(-1.5692697) q[2];
sx q[2];
rz(-2.7806139) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(2.952164) q[1];
sx q[1];
rz(-1.4888114) q[1];
sx q[1];
rz(-1.5174972) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.030121) q[3];
sx q[3];
rz(-2.432593) q[3];
sx q[3];
rz(1.0741155) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-1.063529) q[2];
sx q[2];
rz(-2.0827796) q[2];
sx q[2];
rz(0.25827363) q[2];
rz(0.20049788) q[3];
sx q[3];
rz(-2.343488) q[3];
sx q[3];
rz(2.958278) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9789326) q[0];
sx q[0];
rz(-1.7698092) q[0];
sx q[0];
rz(0.3072511) q[0];
rz(1.9303253) q[1];
sx q[1];
rz(-2.6478421) q[1];
sx q[1];
rz(2.5603851) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2659822) q[0];
sx q[0];
rz(-1.6020244) q[0];
sx q[0];
rz(-1.9755548) q[0];
rz(-pi) q[1];
rz(-2.4318784) q[2];
sx q[2];
rz(-1.6953371) q[2];
sx q[2];
rz(2.1309851) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(2.0705659) q[1];
sx q[1];
rz(-0.80447996) q[1];
sx q[1];
rz(0.2792554) q[1];
rz(-1.2138661) q[3];
sx q[3];
rz(-1.0634616) q[3];
sx q[3];
rz(1.9510253) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-1.4622978) q[2];
sx q[2];
rz(-1.1335979) q[2];
sx q[2];
rz(-3.1104258) q[2];
rz(-2.9160685) q[3];
sx q[3];
rz(-1.3616819) q[3];
sx q[3];
rz(1.0303191) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(pi/2) q[3];
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
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.088260055) q[0];
sx q[0];
rz(-3.0971165) q[0];
sx q[0];
rz(-0.70892507) q[0];
rz(-0.21253474) q[1];
sx q[1];
rz(-2.1268851) q[1];
sx q[1];
rz(-2.2841891) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.014985059) q[0];
sx q[0];
rz(-1.7005655) q[0];
sx q[0];
rz(-3.1364822) q[0];
x q[1];
rz(2.6231758) q[2];
sx q[2];
rz(-0.5974434) q[2];
sx q[2];
rz(-3.0068676) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-1.9136476) q[1];
sx q[1];
rz(-2.1008607) q[1];
sx q[1];
rz(1.1157406) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.69714947) q[3];
sx q[3];
rz(-1.7350041) q[3];
sx q[3];
rz(-0.96443161) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(1.4400441) q[2];
sx q[2];
rz(-2.7893119) q[2];
sx q[2];
rz(2.3360543) q[2];
rz(-2.7696179) q[3];
sx q[3];
rz(-1.493908) q[3];
sx q[3];
rz(2.7443547) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[2];
sx q[3];
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
rz(-3.1086248) q[0];
sx q[0];
rz(-1.2297577) q[0];
sx q[0];
rz(0.97880542) q[0];
rz(-0.72273123) q[1];
sx q[1];
rz(-2.0131854) q[1];
sx q[1];
rz(-0.61789787) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.582983) q[0];
sx q[0];
rz(-1.8454058) q[0];
sx q[0];
rz(-0.78297575) q[0];
rz(-pi) q[1];
rz(2.482867) q[2];
sx q[2];
rz(-1.1470231) q[2];
sx q[2];
rz(1.2032594) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-2.5717575) q[1];
sx q[1];
rz(-1.7392842) q[1];
sx q[1];
rz(1.7525826) q[1];
rz(-0.79094751) q[3];
sx q[3];
rz(-3.0229048) q[3];
sx q[3];
rz(0.044767901) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-1.2021947) q[2];
sx q[2];
rz(-1.3556182) q[2];
sx q[2];
rz(0.00016577684) q[2];
rz(0.58445066) q[3];
sx q[3];
rz(-2.1355459) q[3];
sx q[3];
rz(0.59529006) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7547739) q[0];
sx q[0];
rz(-1.5703572) q[0];
sx q[0];
rz(1.5686709) q[0];
rz(1.8008925) q[1];
sx q[1];
rz(-1.1005713) q[1];
sx q[1];
rz(1.5250199) q[1];
rz(2.1519221) q[2];
sx q[2];
rz(-0.72643092) q[2];
sx q[2];
rz(0.55813172) q[2];
rz(-2.4711193) q[3];
sx q[3];
rz(-2.5994876) q[3];
sx q[3];
rz(2.2710298) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
