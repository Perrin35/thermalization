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
rz(-0.67396069) q[0];
rz(-0.31710467) q[1];
sx q[1];
rz(-1.6333406) q[1];
sx q[1];
rz(-0.79467264) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3622409) q[0];
sx q[0];
rz(-1.9277713) q[0];
sx q[0];
rz(2.7428521) q[0];
x q[1];
rz(2.652466) q[2];
sx q[2];
rz(-1.7009652) q[2];
sx q[2];
rz(-0.39718539) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(1.7197345) q[1];
sx q[1];
rz(-1.0654447) q[1];
sx q[1];
rz(-1.1769597) q[1];
x q[2];
rz(-3.0029293) q[3];
sx q[3];
rz(-1.1274459) q[3];
sx q[3];
rz(-1.5594304) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-0.47544605) q[2];
sx q[2];
rz(-1.1300766) q[2];
sx q[2];
rz(-1.8165992) q[2];
rz(0.31630668) q[3];
sx q[3];
rz(-0.91747147) q[3];
sx q[3];
rz(-0.90707183) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[3];
sx q[3];
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
rz(1.1871724) q[0];
sx q[0];
rz(-0.81048727) q[0];
sx q[0];
rz(0.95970884) q[0];
rz(-1.6540487) q[1];
sx q[1];
rz(-2.4580749) q[1];
sx q[1];
rz(-2.5141292) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.42694416) q[0];
sx q[0];
rz(-1.5776002) q[0];
sx q[0];
rz(3.1338918) q[0];
rz(-pi) q[1];
x q[1];
rz(1.5928028) q[2];
sx q[2];
rz(-1.7823879) q[2];
sx q[2];
rz(1.6197268) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-1.5261425) q[1];
sx q[1];
rz(-2.8612988) q[1];
sx q[1];
rz(-1.0326833) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.5454607) q[3];
sx q[3];
rz(-1.5583894) q[3];
sx q[3];
rz(2.5125463) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-2.6015357) q[2];
sx q[2];
rz(-1.7404218) q[2];
sx q[2];
rz(2.1616948) q[2];
rz(1.4510138) q[3];
sx q[3];
rz(-2.4308949) q[3];
sx q[3];
rz(-1.6987945) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5805805) q[0];
sx q[0];
rz(-0.0479527) q[0];
sx q[0];
rz(-1.3797492) q[0];
rz(-0.17669949) q[1];
sx q[1];
rz(-1.9243762) q[1];
sx q[1];
rz(1.3476936) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.1130226) q[0];
sx q[0];
rz(-1.1060113) q[0];
sx q[0];
rz(-1.3010498) q[0];
rz(-0.10920306) q[2];
sx q[2];
rz(-2.9401527) q[2];
sx q[2];
rz(0.86004721) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(1.2587535) q[1];
sx q[1];
rz(-1.4708609) q[1];
sx q[1];
rz(-0.12211166) q[1];
rz(-pi) q[2];
x q[2];
rz(1.1191145) q[3];
sx q[3];
rz(-0.20483769) q[3];
sx q[3];
rz(-1.1289489) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-1.1220864) q[2];
sx q[2];
rz(-2.2821189) q[2];
sx q[2];
rz(-1.1536095) q[2];
rz(-2.7518318) q[3];
sx q[3];
rz(-0.4699769) q[3];
sx q[3];
rz(-0.78891689) q[3];
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
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6692114) q[0];
sx q[0];
rz(-0.62494576) q[0];
sx q[0];
rz(-0.56458449) q[0];
rz(-3.0584884) q[1];
sx q[1];
rz(-1.233498) q[1];
sx q[1];
rz(2.9262503) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4190061) q[0];
sx q[0];
rz(-2.3424087) q[0];
sx q[0];
rz(-2.0439842) q[0];
rz(-pi) q[1];
x q[1];
rz(0.16673659) q[2];
sx q[2];
rz(-0.98966375) q[2];
sx q[2];
rz(2.2819448) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-1.5147869) q[1];
sx q[1];
rz(-1.6759911) q[1];
sx q[1];
rz(1.3749966) q[1];
x q[2];
rz(-2.4550291) q[3];
sx q[3];
rz(-2.4719704) q[3];
sx q[3];
rz(-1.3470105) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(2.5061491) q[2];
sx q[2];
rz(-1.0887257) q[2];
sx q[2];
rz(-0.60194683) q[2];
rz(2.4510032) q[3];
sx q[3];
rz(-0.75406051) q[3];
sx q[3];
rz(-0.62600342) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
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
rz(-0.85201207) q[0];
sx q[0];
rz(-2.1471725) q[0];
sx q[0];
rz(0.96486282) q[0];
rz(-3.124974) q[1];
sx q[1];
rz(-2.7162572) q[1];
sx q[1];
rz(0.66666493) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.71320888) q[0];
sx q[0];
rz(-0.41233006) q[0];
sx q[0];
rz(1.8250188) q[0];
x q[1];
rz(-2.1140852) q[2];
sx q[2];
rz(-2.5230061) q[2];
sx q[2];
rz(1.0416043) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(2.7409089) q[1];
sx q[1];
rz(-0.92036696) q[1];
sx q[1];
rz(1.9294192) q[1];
x q[2];
rz(-2.9186451) q[3];
sx q[3];
rz(-1.8604391) q[3];
sx q[3];
rz(1.851561) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-2.2255286) q[2];
sx q[2];
rz(-1.3801882) q[2];
sx q[2];
rz(-1.8886245) q[2];
rz(1.8270252) q[3];
sx q[3];
rz(-1.9606749) q[3];
sx q[3];
rz(-1.130828) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
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
rz(2.0261633) q[0];
sx q[0];
rz(-0.9342497) q[0];
sx q[0];
rz(-1.7793659) q[0];
rz(-1.399614) q[1];
sx q[1];
rz(-1.5618961) q[1];
sx q[1];
rz(1.5302352) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0936733) q[0];
sx q[0];
rz(-1.1512655) q[0];
sx q[0];
rz(1.490926) q[0];
x q[1];
rz(-2.5189581) q[2];
sx q[2];
rz(-0.52402516) q[2];
sx q[2];
rz(-0.084409075) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(1.02009) q[1];
sx q[1];
rz(-1.1768747) q[1];
sx q[1];
rz(-2.9083088) q[1];
rz(3.0909782) q[3];
sx q[3];
rz(-1.8236056) q[3];
sx q[3];
rz(1.8265754) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-0.35649148) q[2];
sx q[2];
rz(-1.0078398) q[2];
sx q[2];
rz(3.1121758) q[2];
rz(-1.6507089) q[3];
sx q[3];
rz(-0.95508948) q[3];
sx q[3];
rz(1.9194549) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.71762639) q[0];
sx q[0];
rz(-2.3762149) q[0];
sx q[0];
rz(-2.9587865) q[0];
rz(-0.43552661) q[1];
sx q[1];
rz(-0.85414779) q[1];
sx q[1];
rz(-0.51087728) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2832527) q[0];
sx q[0];
rz(-2.8562299) q[0];
sx q[0];
rz(0.39985379) q[0];
rz(2.4099318) q[2];
sx q[2];
rz(-0.96207843) q[2];
sx q[2];
rz(2.4351956) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.70799815) q[1];
sx q[1];
rz(-0.79213789) q[1];
sx q[1];
rz(2.5670693) q[1];
rz(-pi) q[2];
rz(-2.0931582) q[3];
sx q[3];
rz(-1.8494542) q[3];
sx q[3];
rz(-3.0184961) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(1.9110979) q[2];
sx q[2];
rz(-2.348933) q[2];
sx q[2];
rz(0.62057173) q[2];
rz(-2.7190322) q[3];
sx q[3];
rz(-1.3603323) q[3];
sx q[3];
rz(1.2703936) q[3];
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
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.34185103) q[0];
sx q[0];
rz(-2.7809679) q[0];
sx q[0];
rz(-1.0890695) q[0];
rz(-2.0607121) q[1];
sx q[1];
rz(-1.2842478) q[1];
sx q[1];
rz(-0.63527766) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.23112049) q[0];
sx q[0];
rz(-1.4164682) q[0];
sx q[0];
rz(1.9327823) q[0];
x q[1];
rz(1.073457) q[2];
sx q[2];
rz(-1.724223) q[2];
sx q[2];
rz(1.9286326) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-1.2891846) q[1];
sx q[1];
rz(-1.8361366) q[1];
sx q[1];
rz(-0.62383382) q[1];
rz(0.20181228) q[3];
sx q[3];
rz(-1.0045856) q[3];
sx q[3];
rz(1.7081529) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-0.13715956) q[2];
sx q[2];
rz(-1.4617504) q[2];
sx q[2];
rz(0.00017246406) q[2];
rz(0.6862644) q[3];
sx q[3];
rz(-0.72824794) q[3];
sx q[3];
rz(0.85787684) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8453318) q[0];
sx q[0];
rz(-0.91264549) q[0];
sx q[0];
rz(-2.4849179) q[0];
rz(0.36390057) q[1];
sx q[1];
rz(-0.44635043) q[1];
sx q[1];
rz(2.231853) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1247524) q[0];
sx q[0];
rz(-2.7144055) q[0];
sx q[0];
rz(0.70144002) q[0];
x q[1];
rz(-0.54833834) q[2];
sx q[2];
rz(-0.40204918) q[2];
sx q[2];
rz(2.6532432) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-0.39975702) q[1];
sx q[1];
rz(-1.4474807) q[1];
sx q[1];
rz(2.9305305) q[1];
rz(-pi) q[2];
rz(-2.2496201) q[3];
sx q[3];
rz(-2.744305) q[3];
sx q[3];
rz(0.15115034) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-0.83071128) q[2];
sx q[2];
rz(-2.723366) q[2];
sx q[2];
rz(-0.48661423) q[2];
rz(-3.0330372) q[3];
sx q[3];
rz(-2.4311549) q[3];
sx q[3];
rz(1.2238812) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.20031032) q[0];
sx q[0];
rz(-1.0815942) q[0];
sx q[0];
rz(-2.6249028) q[0];
rz(-1.5453045) q[1];
sx q[1];
rz(-2.089112) q[1];
sx q[1];
rz(-2.2470078) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0588194) q[0];
sx q[0];
rz(-1.528577) q[0];
sx q[0];
rz(-0.93725462) q[0];
rz(2.0403433) q[2];
sx q[2];
rz(-1.38648) q[2];
sx q[2];
rz(1.7061526) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.6539508) q[1];
sx q[1];
rz(-2.7751121) q[1];
sx q[1];
rz(-1.3932863) q[1];
x q[2];
rz(2.9856332) q[3];
sx q[3];
rz(-0.92658639) q[3];
sx q[3];
rz(2.9498266) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-0.14780012) q[2];
sx q[2];
rz(-2.1788969) q[2];
sx q[2];
rz(-2.1195892) q[2];
rz(1.3379124) q[3];
sx q[3];
rz(-1.6513848) q[3];
sx q[3];
rz(-2.5125304) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.548303) q[0];
sx q[0];
rz(-1.5240482) q[0];
sx q[0];
rz(-0.89292009) q[0];
rz(-1.0409566) q[1];
sx q[1];
rz(-3.0491842) q[1];
sx q[1];
rz(-1.4728117) q[1];
rz(-0.9085761) q[2];
sx q[2];
rz(-2.3618345) q[2];
sx q[2];
rz(0.97710412) q[2];
rz(-1.8716807) q[3];
sx q[3];
rz(-1.5242929) q[3];
sx q[3];
rz(2.268189) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
