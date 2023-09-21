OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-0.43006858) q[0];
sx q[0];
rz(-3.0741337) q[0];
sx q[0];
rz(2.467632) q[0];
rz(-0.31710467) q[1];
sx q[1];
rz(-1.6333406) q[1];
sx q[1];
rz(-0.79467264) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0792102) q[0];
sx q[0];
rz(-1.198472) q[0];
sx q[0];
rz(-1.9553493) q[0];
rz(-pi) q[1];
rz(-1.7180213) q[2];
sx q[2];
rz(-2.0554254) q[2];
sx q[2];
rz(2.0369612) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-0.049584576) q[1];
sx q[1];
rz(-1.2283748) q[1];
sx q[1];
rz(2.601806) q[1];
x q[2];
rz(-1.2875597) q[3];
sx q[3];
rz(-0.4631511) q[3];
sx q[3];
rz(1.8739665) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(2.6661466) q[2];
sx q[2];
rz(-1.1300766) q[2];
sx q[2];
rz(1.3249935) q[2];
rz(2.825286) q[3];
sx q[3];
rz(-2.2241212) q[3];
sx q[3];
rz(-0.90707183) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
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
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1871724) q[0];
sx q[0];
rz(-0.81048727) q[0];
sx q[0];
rz(-0.95970884) q[0];
rz(1.487544) q[1];
sx q[1];
rz(-0.68351775) q[1];
sx q[1];
rz(-0.62746343) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7213631) q[0];
sx q[0];
rz(-0.010275928) q[0];
sx q[0];
rz(-2.4179439) q[0];
rz(-pi) q[1];
rz(1.5928028) q[2];
sx q[2];
rz(-1.3592048) q[2];
sx q[2];
rz(1.5218658) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(1.5261425) q[1];
sx q[1];
rz(-0.28029385) q[1];
sx q[1];
rz(-1.0326833) q[1];
x q[2];
rz(2.0262358) q[3];
sx q[3];
rz(-0.028209837) q[3];
sx q[3];
rz(-1.7445604) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(2.6015357) q[2];
sx q[2];
rz(-1.4011708) q[2];
sx q[2];
rz(-0.97989782) q[2];
rz(1.4510138) q[3];
sx q[3];
rz(-2.4308949) q[3];
sx q[3];
rz(1.4427982) q[3];
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
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5805805) q[0];
sx q[0];
rz(-0.0479527) q[0];
sx q[0];
rz(-1.3797492) q[0];
rz(0.17669949) q[1];
sx q[1];
rz(-1.2172164) q[1];
sx q[1];
rz(-1.7938991) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7226669) q[0];
sx q[0];
rz(-1.8113266) q[0];
sx q[0];
rz(0.47970432) q[0];
rz(-pi) q[1];
rz(0.20027192) q[2];
sx q[2];
rz(-1.592604) q[2];
sx q[2];
rz(-2.3238316) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(0.29979953) q[1];
sx q[1];
rz(-1.6922957) q[1];
sx q[1];
rz(-1.4701162) q[1];
x q[2];
rz(-2.0224781) q[3];
sx q[3];
rz(-0.20483769) q[3];
sx q[3];
rz(-1.1289489) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-2.0195062) q[2];
sx q[2];
rz(-2.2821189) q[2];
sx q[2];
rz(-1.9879831) q[2];
rz(-2.7518318) q[3];
sx q[3];
rz(-0.4699769) q[3];
sx q[3];
rz(2.3526758) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.47238123) q[0];
sx q[0];
rz(-2.5166469) q[0];
sx q[0];
rz(-2.5770082) q[0];
rz(3.0584884) q[1];
sx q[1];
rz(-1.233498) q[1];
sx q[1];
rz(0.21534236) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.089245361) q[0];
sx q[0];
rz(-0.87886341) q[0];
sx q[0];
rz(2.7034876) q[0];
rz(-pi) q[1];
x q[1];
rz(2.9748561) q[2];
sx q[2];
rz(-0.98966375) q[2];
sx q[2];
rz(-2.2819448) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-1.5147869) q[1];
sx q[1];
rz(-1.4656015) q[1];
sx q[1];
rz(1.766596) q[1];
rz(-pi) q[2];
rz(-2.5921949) q[3];
sx q[3];
rz(-1.975172) q[3];
sx q[3];
rz(-2.3466563) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-0.6354436) q[2];
sx q[2];
rz(-1.0887257) q[2];
sx q[2];
rz(-0.60194683) q[2];
rz(-0.69058949) q[3];
sx q[3];
rz(-2.3875321) q[3];
sx q[3];
rz(0.62600342) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.85201207) q[0];
sx q[0];
rz(-0.99442011) q[0];
sx q[0];
rz(0.96486282) q[0];
rz(3.124974) q[1];
sx q[1];
rz(-0.42533541) q[1];
sx q[1];
rz(0.66666493) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.71320888) q[0];
sx q[0];
rz(-2.7292626) q[0];
sx q[0];
rz(1.3165738) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.0275074) q[2];
sx q[2];
rz(-0.61858656) q[2];
sx q[2];
rz(-2.0999883) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(2.1866245) q[1];
sx q[1];
rz(-0.72997) q[1];
sx q[1];
rz(-0.43220046) q[1];
rz(-pi) q[2];
rz(-2.9186451) q[3];
sx q[3];
rz(-1.8604391) q[3];
sx q[3];
rz(1.851561) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-0.91606402) q[2];
sx q[2];
rz(-1.3801882) q[2];
sx q[2];
rz(1.8886245) q[2];
rz(-1.3145674) q[3];
sx q[3];
rz(-1.1809177) q[3];
sx q[3];
rz(-2.0107646) q[3];
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
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1154293) q[0];
sx q[0];
rz(-0.9342497) q[0];
sx q[0];
rz(1.7793659) q[0];
rz(1.399614) q[1];
sx q[1];
rz(-1.5618961) q[1];
sx q[1];
rz(1.6113575) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0936733) q[0];
sx q[0];
rz(-1.9903272) q[0];
sx q[0];
rz(-1.490926) q[0];
rz(-0.43892626) q[2];
sx q[2];
rz(-1.274684) q[2];
sx q[2];
rz(2.211328) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(2.6758319) q[1];
sx q[1];
rz(-0.45468802) q[1];
sx q[1];
rz(2.0783706) q[1];
rz(-pi) q[2];
x q[2];
rz(0.050614428) q[3];
sx q[3];
rz(-1.8236056) q[3];
sx q[3];
rz(-1.8265754) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(0.35649148) q[2];
sx q[2];
rz(-1.0078398) q[2];
sx q[2];
rz(0.029416857) q[2];
rz(1.6507089) q[3];
sx q[3];
rz(-2.1865032) q[3];
sx q[3];
rz(1.9194549) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
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
rz(-2.706066) q[1];
sx q[1];
rz(-0.85414779) q[1];
sx q[1];
rz(0.51087728) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.85833997) q[0];
sx q[0];
rz(-2.8562299) q[0];
sx q[0];
rz(0.39985379) q[0];
x q[1];
rz(-2.3235544) q[2];
sx q[2];
rz(-2.15089) q[2];
sx q[2];
rz(1.3387418) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(0.70799815) q[1];
sx q[1];
rz(-0.79213789) q[1];
sx q[1];
rz(-0.57452332) q[1];
rz(-pi) q[2];
rz(1.0501409) q[3];
sx q[3];
rz(-2.5556892) q[3];
sx q[3];
rz(-1.0019765) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-1.2304948) q[2];
sx q[2];
rz(-0.7926597) q[2];
sx q[2];
rz(-0.62057173) q[2];
rz(0.42256045) q[3];
sx q[3];
rz(-1.7812604) q[3];
sx q[3];
rz(1.871199) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
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
rz(-0.34185103) q[0];
sx q[0];
rz(-0.36062476) q[0];
sx q[0];
rz(2.0525232) q[0];
rz(-2.0607121) q[1];
sx q[1];
rz(-1.2842478) q[1];
sx q[1];
rz(-0.63527766) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7437744) q[0];
sx q[0];
rz(-1.9282856) q[0];
sx q[0];
rz(-2.9767569) q[0];
x q[1];
rz(2.0681357) q[2];
sx q[2];
rz(-1.724223) q[2];
sx q[2];
rz(-1.9286326) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.852408) q[1];
sx q[1];
rz(-1.8361366) q[1];
sx q[1];
rz(-2.5177588) q[1];
x q[2];
rz(1.2653207) q[3];
sx q[3];
rz(-2.5442122) q[3];
sx q[3];
rz(1.7978096) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(0.13715956) q[2];
sx q[2];
rz(-1.4617504) q[2];
sx q[2];
rz(3.1414202) q[2];
rz(2.4553283) q[3];
sx q[3];
rz(-2.4133447) q[3];
sx q[3];
rz(0.85787684) q[3];
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
rz(-pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2962608) q[0];
sx q[0];
rz(-2.2289472) q[0];
sx q[0];
rz(-0.6566748) q[0];
rz(0.36390057) q[1];
sx q[1];
rz(-0.44635043) q[1];
sx q[1];
rz(2.231853) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1014935) q[0];
sx q[0];
rz(-1.8414521) q[0];
sx q[0];
rz(-0.33466848) q[0];
rz(2.5932543) q[2];
sx q[2];
rz(-2.7395435) q[2];
sx q[2];
rz(0.48834947) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-1.9442056) q[1];
sx q[1];
rz(-1.3613609) q[1];
sx q[1];
rz(1.6968813) q[1];
x q[2];
rz(-1.2551366) q[3];
sx q[3];
rz(-1.3254032) q[3];
sx q[3];
rz(1.0823702) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-0.83071128) q[2];
sx q[2];
rz(-0.41822663) q[2];
sx q[2];
rz(-2.6549784) q[2];
rz(0.1085554) q[3];
sx q[3];
rz(-2.4311549) q[3];
sx q[3];
rz(1.2238812) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
x q[2];
rz(pi/2) q[3];
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
rz(2.9412823) q[0];
sx q[0];
rz(-1.0815942) q[0];
sx q[0];
rz(-2.6249028) q[0];
rz(-1.5453045) q[1];
sx q[1];
rz(-2.089112) q[1];
sx q[1];
rz(-2.2470078) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.54297011) q[0];
sx q[0];
rz(-0.93790903) q[0];
sx q[0];
rz(3.0892239) q[0];
x q[1];
rz(1.1012494) q[2];
sx q[2];
rz(-1.7551127) q[2];
sx q[2];
rz(1.7061526) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(2.6539508) q[1];
sx q[1];
rz(-0.36648053) q[1];
sx q[1];
rz(-1.3932863) q[1];
rz(-1.3668725) q[3];
sx q[3];
rz(-2.4813934) q[3];
sx q[3];
rz(3.077293) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(2.9937925) q[2];
sx q[2];
rz(-2.1788969) q[2];
sx q[2];
rz(1.0220035) q[2];
rz(-1.3379124) q[3];
sx q[3];
rz(-1.4902078) q[3];
sx q[3];
rz(-2.5125304) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.548303) q[0];
sx q[0];
rz(-1.6175445) q[0];
sx q[0];
rz(2.2486726) q[0];
rz(2.100636) q[1];
sx q[1];
rz(-3.0491842) q[1];
sx q[1];
rz(-1.4728117) q[1];
rz(-2.5953318) q[2];
sx q[2];
rz(-0.98304521) q[2];
sx q[2];
rz(1.8084768) q[2];
rz(0.048687497) q[3];
sx q[3];
rz(-1.2702474) q[3];
sx q[3];
rz(-2.4297759) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];