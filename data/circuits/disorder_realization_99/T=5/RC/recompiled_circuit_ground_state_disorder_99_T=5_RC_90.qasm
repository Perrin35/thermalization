OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(0.91166624) q[0];
sx q[0];
rz(-0.32719964) q[0];
sx q[0];
rz(1.3258452) q[0];
rz(1.3658547) q[1];
sx q[1];
rz(-0.39443016) q[1];
sx q[1];
rz(0.78864133) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.45730725) q[0];
sx q[0];
rz(-0.72170242) q[0];
sx q[0];
rz(-1.0941605) q[0];
rz(-pi) q[1];
rz(-2.4633258) q[2];
sx q[2];
rz(-2.1951198) q[2];
sx q[2];
rz(-2.2140142) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.44581283) q[1];
sx q[1];
rz(-1.7632293) q[1];
sx q[1];
rz(-2.8532527) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.9427845) q[3];
sx q[3];
rz(-1.4664141) q[3];
sx q[3];
rz(-1.2918772) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-2.6476562) q[2];
sx q[2];
rz(-1.6948573) q[2];
sx q[2];
rz(2.9392865) q[2];
rz(-0.10937396) q[3];
sx q[3];
rz(-0.60905639) q[3];
sx q[3];
rz(-3.0561395) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0626471) q[0];
sx q[0];
rz(-2.2181692) q[0];
sx q[0];
rz(1.9664636) q[0];
rz(-1.6773978) q[1];
sx q[1];
rz(-1.0025832) q[1];
sx q[1];
rz(-0.35951231) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.26665154) q[0];
sx q[0];
rz(-1.4066937) q[0];
sx q[0];
rz(-3.0709188) q[0];
x q[1];
rz(1.9486675) q[2];
sx q[2];
rz(-1.0456351) q[2];
sx q[2];
rz(1.988387) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(1.2613758) q[1];
sx q[1];
rz(-1.2296074) q[1];
sx q[1];
rz(2.5009194) q[1];
rz(-0.80163281) q[3];
sx q[3];
rz(-0.75770658) q[3];
sx q[3];
rz(-3.0996292) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.8466865) q[2];
sx q[2];
rz(-2.0945695) q[2];
sx q[2];
rz(0.27493757) q[2];
rz(-1.8168195) q[3];
sx q[3];
rz(-1.6144269) q[3];
sx q[3];
rz(1.2980609) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0777968) q[0];
sx q[0];
rz(-2.9834788) q[0];
sx q[0];
rz(2.7445444) q[0];
rz(2.5322757) q[1];
sx q[1];
rz(-1.3343697) q[1];
sx q[1];
rz(-1.5416001) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1572874) q[0];
sx q[0];
rz(-1.8616891) q[0];
sx q[0];
rz(0.10910927) q[0];
x q[1];
rz(-3.1290319) q[2];
sx q[2];
rz(-1.505559) q[2];
sx q[2];
rz(3.1029004) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(1.6172495) q[1];
sx q[1];
rz(-1.5876974) q[1];
sx q[1];
rz(0.84417697) q[1];
rz(2.3261689) q[3];
sx q[3];
rz(-2.2924726) q[3];
sx q[3];
rz(-1.4072925) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-2.2043173) q[2];
sx q[2];
rz(-1.5466362) q[2];
sx q[2];
rz(0.23359648) q[2];
rz(-1.8448081) q[3];
sx q[3];
rz(-1.1917043) q[3];
sx q[3];
rz(-2.4209723) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.63491708) q[0];
sx q[0];
rz(-1.1577865) q[0];
sx q[0];
rz(-1.7764212) q[0];
rz(-1.8695976) q[1];
sx q[1];
rz(-1.1993661) q[1];
sx q[1];
rz(-1.8796399) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.015975706) q[0];
sx q[0];
rz(-1.7628094) q[0];
sx q[0];
rz(-0.58108347) q[0];
rz(-pi) q[1];
rz(1.8795525) q[2];
sx q[2];
rz(-1.595904) q[2];
sx q[2];
rz(-1.0756191) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-3.1168667) q[1];
sx q[1];
rz(-2.9477009) q[1];
sx q[1];
rz(2.7996382) q[1];
x q[2];
rz(-0.19660321) q[3];
sx q[3];
rz(-2.3075309) q[3];
sx q[3];
rz(1.3722668) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.9103553) q[2];
sx q[2];
rz(-2.1114025) q[2];
sx q[2];
rz(-0.78872284) q[2];
rz(-1.2809058) q[3];
sx q[3];
rz(-1.3642045) q[3];
sx q[3];
rz(2.1196712) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7541517) q[0];
sx q[0];
rz(-1.0197637) q[0];
sx q[0];
rz(0.66106853) q[0];
rz(-1.1152274) q[1];
sx q[1];
rz(-1.8293569) q[1];
sx q[1];
rz(2.1818395) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.57058731) q[0];
sx q[0];
rz(-0.59413213) q[0];
sx q[0];
rz(0.49555619) q[0];
rz(2.9767545) q[2];
sx q[2];
rz(-1.0292813) q[2];
sx q[2];
rz(-0.37357003) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-2.5670523) q[1];
sx q[1];
rz(-0.31992074) q[1];
sx q[1];
rz(1.1362651) q[1];
x q[2];
rz(-1.8314535) q[3];
sx q[3];
rz(-2.4803526) q[3];
sx q[3];
rz(-2.6186297) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(0.58012086) q[2];
sx q[2];
rz(-0.46130195) q[2];
sx q[2];
rz(1.0599773) q[2];
rz(-0.75016841) q[3];
sx q[3];
rz(-1.3786022) q[3];
sx q[3];
rz(-1.9157971) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.2724514) q[0];
sx q[0];
rz(-1.5831818) q[0];
sx q[0];
rz(-2.9172752) q[0];
rz(-1.51651) q[1];
sx q[1];
rz(-2.5254011) q[1];
sx q[1];
rz(0.075627653) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6859378) q[0];
sx q[0];
rz(-1.9561531) q[0];
sx q[0];
rz(-2.1598201) q[0];
rz(-pi) q[1];
x q[1];
rz(0.096616726) q[2];
sx q[2];
rz(-1.8732837) q[2];
sx q[2];
rz(-2.5563478) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-0.26201263) q[1];
sx q[1];
rz(-2.0063489) q[1];
sx q[1];
rz(2.0836989) q[1];
rz(-0.43419713) q[3];
sx q[3];
rz(-0.89354529) q[3];
sx q[3];
rz(-1.7627258) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-1.6077967) q[2];
sx q[2];
rz(-2.3888612) q[2];
sx q[2];
rz(0.047018615) q[2];
rz(-2.4646344) q[3];
sx q[3];
rz(-1.8432063) q[3];
sx q[3];
rz(0.025378749) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.15636477) q[0];
sx q[0];
rz(-1.4893091) q[0];
sx q[0];
rz(-1.697668) q[0];
rz(0.9230744) q[1];
sx q[1];
rz(-2.1363027) q[1];
sx q[1];
rz(-0.18327555) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3472206) q[0];
sx q[0];
rz(-2.1275015) q[0];
sx q[0];
rz(1.8666202) q[0];
rz(-pi) q[1];
x q[1];
rz(1.9940111) q[2];
sx q[2];
rz(-1.9789038) q[2];
sx q[2];
rz(1.0427047) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(2.7828655) q[1];
sx q[1];
rz(-1.0762097) q[1];
sx q[1];
rz(-2.7703157) q[1];
x q[2];
rz(-1.9952946) q[3];
sx q[3];
rz(-2.5190926) q[3];
sx q[3];
rz(1.5869035) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.034059374) q[2];
sx q[2];
rz(-1.6851765) q[2];
sx q[2];
rz(0.019406645) q[2];
rz(2.9803045) q[3];
sx q[3];
rz(-2.1262157) q[3];
sx q[3];
rz(1.2279145) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1041383) q[0];
sx q[0];
rz(-2.7327974) q[0];
sx q[0];
rz(0.092305146) q[0];
rz(-1.1760938) q[1];
sx q[1];
rz(-0.25830019) q[1];
sx q[1];
rz(-1.4091122) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4482142) q[0];
sx q[0];
rz(-1.6758159) q[0];
sx q[0];
rz(3.0870221) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.1827132) q[2];
sx q[2];
rz(-1.1772708) q[2];
sx q[2];
rz(-2.3589691) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(0.17615023) q[1];
sx q[1];
rz(-1.9561844) q[1];
sx q[1];
rz(1.8606645) q[1];
rz(-pi) q[2];
rz(-0.73958379) q[3];
sx q[3];
rz(-0.7639262) q[3];
sx q[3];
rz(0.52810625) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.041542355) q[2];
sx q[2];
rz(-1.5771733) q[2];
sx q[2];
rz(1.1735631) q[2];
rz(-0.38343492) q[3];
sx q[3];
rz(-1.3078559) q[3];
sx q[3];
rz(1.0907382) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1293056) q[0];
sx q[0];
rz(-1.6467935) q[0];
sx q[0];
rz(-2.9360085) q[0];
rz(2.3221817) q[1];
sx q[1];
rz(-1.9267547) q[1];
sx q[1];
rz(-0.49088556) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4644829) q[0];
sx q[0];
rz(-2.3620531) q[0];
sx q[0];
rz(-2.5960514) q[0];
x q[1];
rz(-2.728168) q[2];
sx q[2];
rz(-2.4707762) q[2];
sx q[2];
rz(-0.49395032) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(1.6632884) q[1];
sx q[1];
rz(-0.57034796) q[1];
sx q[1];
rz(-1.4823518) q[1];
x q[2];
rz(1.7455467) q[3];
sx q[3];
rz(-0.51769231) q[3];
sx q[3];
rz(-2.2153436) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.72863355) q[2];
sx q[2];
rz(-2.1042991) q[2];
sx q[2];
rz(-1.2566077) q[2];
rz(-2.7058153) q[3];
sx q[3];
rz(-2.0332917) q[3];
sx q[3];
rz(2.2573439) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7213781) q[0];
sx q[0];
rz(-2.2342137) q[0];
sx q[0];
rz(-2.9086928) q[0];
rz(0.13889343) q[1];
sx q[1];
rz(-1.9878309) q[1];
sx q[1];
rz(2.6331666) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.31688894) q[0];
sx q[0];
rz(-2.1068067) q[0];
sx q[0];
rz(-0.30110036) q[0];
x q[1];
rz(2.0783689) q[2];
sx q[2];
rz(-1.8120349) q[2];
sx q[2];
rz(-0.2080179) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-0.99791807) q[1];
sx q[1];
rz(-1.4368205) q[1];
sx q[1];
rz(0.0079946144) q[1];
rz(-pi) q[2];
x q[2];
rz(0.53896972) q[3];
sx q[3];
rz(-0.94851714) q[3];
sx q[3];
rz(-1.301018) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.7306708) q[2];
sx q[2];
rz(-1.6795009) q[2];
sx q[2];
rz(-0.81653583) q[2];
rz(0.38980347) q[3];
sx q[3];
rz(-1.0023508) q[3];
sx q[3];
rz(1.3780814) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.447406) q[0];
sx q[0];
rz(-2.2185855) q[0];
sx q[0];
rz(2.9119281) q[0];
rz(1.5028839) q[1];
sx q[1];
rz(-1.8062183) q[1];
sx q[1];
rz(0.91581215) q[1];
rz(-1.7433132) q[2];
sx q[2];
rz(-1.948922) q[2];
sx q[2];
rz(2.7073635) q[2];
rz(2.9749246) q[3];
sx q[3];
rz(-0.87571908) q[3];
sx q[3];
rz(-2.2805211) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
