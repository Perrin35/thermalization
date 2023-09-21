OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-3.5716612) q[0];
sx q[0];
rz(-0.067458955) q[0];
sx q[0];
rz(10.098739) q[0];
rz(2.824488) q[1];
sx q[1];
rz(-1.5082521) q[1];
sx q[1];
rz(-2.34692) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.062382467) q[0];
sx q[0];
rz(-1.9431207) q[0];
sx q[0];
rz(-1.1862434) q[0];
rz(1.7180213) q[2];
sx q[2];
rz(-1.0861673) q[2];
sx q[2];
rz(-1.1046315) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-1.0104048) q[1];
sx q[1];
rz(-0.63001761) q[1];
sx q[1];
rz(0.6063993) q[1];
rz(-1.2875597) q[3];
sx q[3];
rz(-2.6784416) q[3];
sx q[3];
rz(1.2676261) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-0.47544605) q[2];
sx q[2];
rz(-2.011516) q[2];
sx q[2];
rz(1.8165992) q[2];
rz(-2.825286) q[3];
sx q[3];
rz(-2.2241212) q[3];
sx q[3];
rz(0.90707183) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9544202) q[0];
sx q[0];
rz(-2.3311054) q[0];
sx q[0];
rz(-2.1818838) q[0];
rz(1.487544) q[1];
sx q[1];
rz(-0.68351775) q[1];
sx q[1];
rz(2.5141292) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7213631) q[0];
sx q[0];
rz(-0.010275928) q[0];
sx q[0];
rz(2.4179439) q[0];
rz(-pi) q[1];
rz(1.5928028) q[2];
sx q[2];
rz(-1.3592048) q[2];
sx q[2];
rz(-1.6197268) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(0.97034772) q[1];
sx q[1];
rz(-1.8106318) q[1];
sx q[1];
rz(-2.99511) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.596132) q[3];
sx q[3];
rz(-1.5583894) q[3];
sx q[3];
rz(0.62904639) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.54005694) q[2];
sx q[2];
rz(-1.4011708) q[2];
sx q[2];
rz(0.97989782) q[2];
rz(1.6905789) q[3];
sx q[3];
rz(-2.4308949) q[3];
sx q[3];
rz(1.6987945) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5805805) q[0];
sx q[0];
rz(-0.0479527) q[0];
sx q[0];
rz(-1.7618435) q[0];
rz(0.17669949) q[1];
sx q[1];
rz(-1.9243762) q[1];
sx q[1];
rz(-1.3476936) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5603148) q[0];
sx q[0];
rz(-2.609195) q[0];
sx q[0];
rz(-2.6530932) q[0];
rz(-pi) q[1];
rz(0.10920306) q[2];
sx q[2];
rz(-0.20143992) q[2];
sx q[2];
rz(0.86004721) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(0.29979953) q[1];
sx q[1];
rz(-1.449297) q[1];
sx q[1];
rz(1.4701162) q[1];
rz(-pi) q[2];
x q[2];
rz(-3.0511608) q[3];
sx q[3];
rz(-1.3867497) q[3];
sx q[3];
rz(-0.66891608) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-1.1220864) q[2];
sx q[2];
rz(-0.85947376) q[2];
sx q[2];
rz(-1.9879831) q[2];
rz(0.38976088) q[3];
sx q[3];
rz(-2.6716158) q[3];
sx q[3];
rz(-2.3526758) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[2];
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
rz(0.47238123) q[0];
sx q[0];
rz(-2.5166469) q[0];
sx q[0];
rz(-0.56458449) q[0];
rz(-0.083104221) q[1];
sx q[1];
rz(-1.233498) q[1];
sx q[1];
rz(0.21534236) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9504844) q[0];
sx q[0];
rz(-1.9035625) q[0];
sx q[0];
rz(-0.82975181) q[0];
rz(-pi) q[1];
rz(-1.3232857) q[2];
sx q[2];
rz(-2.5396721) q[2];
sx q[2];
rz(1.9844696) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(0.54289651) q[1];
sx q[1];
rz(-0.22194949) q[1];
sx q[1];
rz(2.0680244) q[1];
x q[2];
rz(0.68656355) q[3];
sx q[3];
rz(-2.4719704) q[3];
sx q[3];
rz(-1.3470105) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-0.6354436) q[2];
sx q[2];
rz(-1.0887257) q[2];
sx q[2];
rz(2.5396458) q[2];
rz(-0.69058949) q[3];
sx q[3];
rz(-0.75406051) q[3];
sx q[3];
rz(2.5155892) q[3];
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
x q[0];
rz(-pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2895806) q[0];
sx q[0];
rz(-0.99442011) q[0];
sx q[0];
rz(-0.96486282) q[0];
rz(-3.124974) q[1];
sx q[1];
rz(-2.7162572) q[1];
sx q[1];
rz(-2.4749277) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.704741) q[0];
sx q[0];
rz(-1.9691103) q[0];
sx q[0];
rz(3.0320291) q[0];
x q[1];
rz(2.1180192) q[2];
sx q[2];
rz(-1.2663411) q[2];
sx q[2];
rz(-0.98642245) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.40068377) q[1];
sx q[1];
rz(-0.92036696) q[1];
sx q[1];
rz(-1.2121735) q[1];
rz(-1.8673709) q[3];
sx q[3];
rz(-1.7843102) q[3];
sx q[3];
rz(-2.7961658) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(2.2255286) q[2];
sx q[2];
rz(-1.7614044) q[2];
sx q[2];
rz(1.2529681) q[2];
rz(-1.3145674) q[3];
sx q[3];
rz(-1.1809177) q[3];
sx q[3];
rz(1.130828) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
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
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0261633) q[0];
sx q[0];
rz(-2.207343) q[0];
sx q[0];
rz(-1.3622267) q[0];
rz(1.399614) q[1];
sx q[1];
rz(-1.5618961) q[1];
sx q[1];
rz(-1.5302352) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8996457) q[0];
sx q[0];
rz(-0.42662222) q[0];
sx q[0];
rz(0.17701478) q[0];
x q[1];
rz(-1.2457232) q[2];
sx q[2];
rz(-1.1522066) q[2];
sx q[2];
rz(-2.3649154) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.64165243) q[1];
sx q[1];
rz(-1.3556726) q[1];
sx q[1];
rz(-1.9745449) q[1];
rz(-1.3176765) q[3];
sx q[3];
rz(-1.6198006) q[3];
sx q[3];
rz(-2.8984836) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(0.35649148) q[2];
sx q[2];
rz(-1.0078398) q[2];
sx q[2];
rz(-3.1121758) q[2];
rz(-1.6507089) q[3];
sx q[3];
rz(-0.95508948) q[3];
sx q[3];
rz(1.9194549) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.71762639) q[0];
sx q[0];
rz(-0.76537776) q[0];
sx q[0];
rz(-0.18280612) q[0];
rz(2.706066) q[1];
sx q[1];
rz(-2.2874449) q[1];
sx q[1];
rz(0.51087728) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.32719192) q[0];
sx q[0];
rz(-1.4609903) q[0];
sx q[0];
rz(-0.26392428) q[0];
x q[1];
rz(2.3350231) q[2];
sx q[2];
rz(-2.2273846) q[2];
sx q[2];
rz(2.8442531) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(2.4335945) q[1];
sx q[1];
rz(-2.3494548) q[1];
sx q[1];
rz(2.5670693) q[1];
rz(1.0501409) q[3];
sx q[3];
rz(-0.58590349) q[3];
sx q[3];
rz(-2.1396162) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(1.9110979) q[2];
sx q[2];
rz(-0.7926597) q[2];
sx q[2];
rz(2.5210209) q[2];
rz(2.7190322) q[3];
sx q[3];
rz(-1.7812604) q[3];
sx q[3];
rz(-1.871199) q[3];
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
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7997416) q[0];
sx q[0];
rz(-0.36062476) q[0];
sx q[0];
rz(1.0890695) q[0];
rz(-2.0607121) q[1];
sx q[1];
rz(-1.2842478) q[1];
sx q[1];
rz(-0.63527766) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1875293) q[0];
sx q[0];
rz(-2.7494193) q[0];
sx q[0];
rz(-1.1568882) q[0];
rz(-1.8842472) q[2];
sx q[2];
rz(-2.6230276) q[2];
sx q[2];
rz(-2.5093362) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(3.0465225) q[1];
sx q[1];
rz(-0.97192837) q[1];
sx q[1];
rz(-1.2477161) q[1];
x q[2];
rz(-0.99526309) q[3];
sx q[3];
rz(-1.4008153) q[3];
sx q[3];
rz(-0.02804027) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-0.13715956) q[2];
sx q[2];
rz(-1.4617504) q[2];
sx q[2];
rz(-3.1414202) q[2];
rz(-0.6862644) q[3];
sx q[3];
rz(-0.72824794) q[3];
sx q[3];
rz(-0.85787684) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2962608) q[0];
sx q[0];
rz(-2.2289472) q[0];
sx q[0];
rz(0.6566748) q[0];
rz(0.36390057) q[1];
sx q[1];
rz(-0.44635043) q[1];
sx q[1];
rz(2.231853) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1014935) q[0];
sx q[0];
rz(-1.8414521) q[0];
sx q[0];
rz(-2.8069242) q[0];
rz(-pi) q[1];
rz(2.7934974) q[2];
sx q[2];
rz(-1.7762134) q[2];
sx q[2];
rz(1.5470488) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(1.197387) q[1];
sx q[1];
rz(-1.7802317) q[1];
sx q[1];
rz(1.4447114) q[1];
x q[2];
rz(-0.89197253) q[3];
sx q[3];
rz(-0.39728764) q[3];
sx q[3];
rz(-2.9904423) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.3108814) q[2];
sx q[2];
rz(-0.41822663) q[2];
sx q[2];
rz(0.48661423) q[2];
rz(0.1085554) q[3];
sx q[3];
rz(-0.71043772) q[3];
sx q[3];
rz(-1.2238812) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.20031032) q[0];
sx q[0];
rz(-1.0815942) q[0];
sx q[0];
rz(0.51668984) q[0];
rz(-1.5962881) q[1];
sx q[1];
rz(-1.0524806) q[1];
sx q[1];
rz(0.89458481) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0588194) q[0];
sx q[0];
rz(-1.6130157) q[0];
sx q[0];
rz(-2.204338) q[0];
rz(-1.1012494) q[2];
sx q[2];
rz(-1.7551127) q[2];
sx q[2];
rz(-1.7061526) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-0.48764187) q[1];
sx q[1];
rz(-0.36648053) q[1];
sx q[1];
rz(-1.3932863) q[1];
x q[2];
rz(1.7747202) q[3];
sx q[3];
rz(-2.4813934) q[3];
sx q[3];
rz(3.077293) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(2.9937925) q[2];
sx q[2];
rz(-2.1788969) q[2];
sx q[2];
rz(-2.1195892) q[2];
rz(-1.3379124) q[3];
sx q[3];
rz(-1.4902078) q[3];
sx q[3];
rz(0.62906229) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5932896) q[0];
sx q[0];
rz(-1.5240482) q[0];
sx q[0];
rz(-0.89292009) q[0];
rz(-2.100636) q[1];
sx q[1];
rz(-0.092408471) q[1];
sx q[1];
rz(1.668781) q[1];
rz(0.54626089) q[2];
sx q[2];
rz(-0.98304521) q[2];
sx q[2];
rz(1.8084768) q[2];
rz(1.7265504) q[3];
sx q[3];
rz(-2.8372436) q[3];
sx q[3];
rz(-2.5929034) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
