OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-0.92276031) q[0];
sx q[0];
rz(3.2050686) q[0];
sx q[0];
rz(10.795074) q[0];
rz(0.81075794) q[1];
sx q[1];
rz(-1.4501362) q[1];
sx q[1];
rz(0.22051799) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6761326) q[0];
sx q[0];
rz(-1.6384083) q[0];
sx q[0];
rz(-2.3071482) q[0];
rz(2.8120997) q[2];
sx q[2];
rz(-1.1348083) q[2];
sx q[2];
rz(1.5583676) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(1.619239) q[1];
sx q[1];
rz(-1.5597889) q[1];
sx q[1];
rz(0.22214823) q[1];
rz(-pi) q[2];
rz(2.9412782) q[3];
sx q[3];
rz(-1.2216) q[3];
sx q[3];
rz(2.0053476) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-0.084879547) q[2];
sx q[2];
rz(-3.1352545) q[2];
sx q[2];
rz(-0.81981266) q[2];
rz(0.020717185) q[3];
sx q[3];
rz(-0.31532225) q[3];
sx q[3];
rz(1.0516385) q[3];
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
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.76899511) q[0];
sx q[0];
rz(-0.1854493) q[0];
sx q[0];
rz(0.73082596) q[0];
rz(1.7164879) q[1];
sx q[1];
rz(-2.4162879) q[1];
sx q[1];
rz(-2.6740429) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0529405) q[0];
sx q[0];
rz(-3.0346617) q[0];
sx q[0];
rz(-1.2783952) q[0];
rz(-pi) q[1];
rz(-0.40981648) q[2];
sx q[2];
rz(-1.6011607) q[2];
sx q[2];
rz(1.8259468) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(1.4380355) q[1];
sx q[1];
rz(-0.86895567) q[1];
sx q[1];
rz(-1.7662394) q[1];
x q[2];
rz(2.3043724) q[3];
sx q[3];
rz(-1.3488111) q[3];
sx q[3];
rz(2.8573621) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(1.1656437) q[2];
sx q[2];
rz(-0.011198137) q[2];
sx q[2];
rz(-2.8051918) q[2];
rz(-2.8339556) q[3];
sx q[3];
rz(-0.010713723) q[3];
sx q[3];
rz(-0.78905869) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.12863185) q[0];
sx q[0];
rz(-0.19879453) q[0];
sx q[0];
rz(3.0096753) q[0];
rz(-1.7031274) q[1];
sx q[1];
rz(-1.4142282) q[1];
sx q[1];
rz(1.6579113) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0795201) q[0];
sx q[0];
rz(-1.8490845) q[0];
sx q[0];
rz(2.5132283) q[0];
x q[1];
rz(2.5222579) q[2];
sx q[2];
rz(-3.1115981) q[2];
sx q[2];
rz(-0.66823792) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-3.0358932) q[1];
sx q[1];
rz(-1.4669682) q[1];
sx q[1];
rz(0.052171252) q[1];
rz(-pi) q[2];
rz(-3.0555757) q[3];
sx q[3];
rz(-1.6413851) q[3];
sx q[3];
rz(3.1215145) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.57054532) q[2];
sx q[2];
rz(-1.3344301) q[2];
sx q[2];
rz(1.4578311) q[2];
rz(2.3633862) q[3];
sx q[3];
rz(-0.0056548803) q[3];
sx q[3];
rz(-2.9290504) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(pi/2) q[3];
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
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1218162) q[0];
sx q[0];
rz(-1.7541405) q[0];
sx q[0];
rz(-0.29086599) q[0];
rz(1.5930755) q[1];
sx q[1];
rz(-2.3389108) q[1];
sx q[1];
rz(-3.1150418) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0850459) q[0];
sx q[0];
rz(-1.9246058) q[0];
sx q[0];
rz(-0.36713552) q[0];
rz(-pi) q[1];
rz(3.0316684) q[2];
sx q[2];
rz(-0.4532686) q[2];
sx q[2];
rz(2.088415) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(1.2902879) q[1];
sx q[1];
rz(-2.0495702) q[1];
sx q[1];
rz(-2.4393875) q[1];
x q[2];
rz(-1.7351355) q[3];
sx q[3];
rz(-3.1222635) q[3];
sx q[3];
rz(-0.011474495) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.87963858) q[2];
sx q[2];
rz(-0.018435437) q[2];
sx q[2];
rz(-0.43799841) q[2];
rz(-2.0336464) q[3];
sx q[3];
rz(-3.1236533) q[3];
sx q[3];
rz(1.3560791) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.22439013) q[0];
sx q[0];
rz(-2.6438535) q[0];
sx q[0];
rz(0.14601953) q[0];
rz(3.0085425) q[1];
sx q[1];
rz(-2.0818043) q[1];
sx q[1];
rz(-0.45566794) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3489482) q[0];
sx q[0];
rz(-1.6647208) q[0];
sx q[0];
rz(1.5507038) q[0];
rz(-pi) q[1];
x q[1];
rz(1.8784166) q[2];
sx q[2];
rz(-2.1407749) q[2];
sx q[2];
rz(2.7893699) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-0.6745488) q[1];
sx q[1];
rz(-1.70948) q[1];
sx q[1];
rz(-1.8013181) q[1];
rz(-pi) q[2];
x q[2];
rz(3.0751441) q[3];
sx q[3];
rz(-0.22804582) q[3];
sx q[3];
rz(-1.4888637) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(1.8749775) q[2];
sx q[2];
rz(-3.1217323) q[2];
sx q[2];
rz(-1.88545) q[2];
rz(-2.6202294) q[3];
sx q[3];
rz(-0.25786906) q[3];
sx q[3];
rz(-0.34463394) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1954023) q[0];
sx q[0];
rz(-2.7552216) q[0];
sx q[0];
rz(-0.56889164) q[0];
rz(2.9733114) q[1];
sx q[1];
rz(-1.556309) q[1];
sx q[1];
rz(3.1313484) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5880347) q[0];
sx q[0];
rz(-1.7419683) q[0];
sx q[0];
rz(-1.6200911) q[0];
rz(-pi) q[1];
rz(3.1222759) q[2];
sx q[2];
rz(-1.7512808) q[2];
sx q[2];
rz(-1.3767124) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(3.0371984) q[1];
sx q[1];
rz(-1.5362829) q[1];
sx q[1];
rz(-1.6134279) q[1];
x q[2];
rz(-1.5404352) q[3];
sx q[3];
rz(-0.93685461) q[3];
sx q[3];
rz(-2.674814) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(2.9804618) q[2];
sx q[2];
rz(-2.9176596) q[2];
sx q[2];
rz(2.0437608) q[2];
rz(0.3499507) q[3];
sx q[3];
rz(-1.8643458) q[3];
sx q[3];
rz(2.2866975) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2646645) q[0];
sx q[0];
rz(-2.9783037) q[0];
sx q[0];
rz(0.40633416) q[0];
rz(1.3111275) q[1];
sx q[1];
rz(-2.6268112) q[1];
sx q[1];
rz(-3.0313671) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.742934) q[0];
sx q[0];
rz(-3.0386819) q[0];
sx q[0];
rz(-0.2700996) q[0];
rz(-pi) q[1];
rz(1.5761779) q[2];
sx q[2];
rz(-1.5670243) q[2];
sx q[2];
rz(1.9413858) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-0.87482968) q[1];
sx q[1];
rz(-1.5620666) q[1];
sx q[1];
rz(-1.1432024) q[1];
rz(-pi) q[2];
rz(2.4065287) q[3];
sx q[3];
rz(-2.1218315) q[3];
sx q[3];
rz(0.50526528) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-1.8608287) q[2];
sx q[2];
rz(-0.0088366652) q[2];
sx q[2];
rz(-0.78738085) q[2];
rz(2.8619838) q[3];
sx q[3];
rz(-0.044737261) q[3];
sx q[3];
rz(0.71981049) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(-pi/2) q[3];
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
rz(-0.0019919458) q[0];
sx q[0];
rz(-2.3878492) q[0];
sx q[0];
rz(-1.4308223) q[0];
rz(-0.17499533) q[1];
sx q[1];
rz(-0.7376968) q[1];
sx q[1];
rz(-1.7134679) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5663877) q[0];
sx q[0];
rz(-3.1194127) q[0];
sx q[0];
rz(1.7087858) q[0];
x q[1];
rz(1.578184) q[2];
sx q[2];
rz(-1.5748757) q[2];
sx q[2];
rz(-2.5830327) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-0.54970804) q[1];
sx q[1];
rz(-1.8401405) q[1];
sx q[1];
rz(2.9742105) q[1];
rz(-pi) q[2];
x q[2];
rz(0.28167991) q[3];
sx q[3];
rz(-1.7865281) q[3];
sx q[3];
rz(0.40116596) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-2.3565107) q[2];
sx q[2];
rz(-2.3729237) q[2];
sx q[2];
rz(-1.53995) q[2];
rz(0.29940638) q[3];
sx q[3];
rz(-1.6573903) q[3];
sx q[3];
rz(0.33777344) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.47141075) q[0];
sx q[0];
rz(-3.1306559) q[0];
sx q[0];
rz(2.6543044) q[0];
rz(1.4393073) q[1];
sx q[1];
rz(-0.7936365) q[1];
sx q[1];
rz(-2.7557441) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.841087) q[0];
sx q[0];
rz(-0.39751507) q[0];
sx q[0];
rz(0.15480283) q[0];
x q[1];
rz(-3.0294777) q[2];
sx q[2];
rz(-1.1476074) q[2];
sx q[2];
rz(-1.2982757) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(0.072695168) q[1];
sx q[1];
rz(-2.4494315) q[1];
sx q[1];
rz(3.0897184) q[1];
rz(-pi) q[2];
rz(-0.33521426) q[3];
sx q[3];
rz(-1.5371369) q[3];
sx q[3];
rz(1.0531283) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(1.7172598) q[2];
sx q[2];
rz(-3.11185) q[2];
sx q[2];
rz(2.7359803) q[2];
rz(2.3154216) q[3];
sx q[3];
rz(-3.0982389) q[3];
sx q[3];
rz(2.1521173) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9721603) q[0];
sx q[0];
rz(-2.5763474) q[0];
sx q[0];
rz(-1.3954847) q[0];
rz(-1.6204429) q[1];
sx q[1];
rz(-0.65137678) q[1];
sx q[1];
rz(1.3491389) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5573709) q[0];
sx q[0];
rz(-0.85600805) q[0];
sx q[0];
rz(-2.8609311) q[0];
x q[1];
rz(-0.17041309) q[2];
sx q[2];
rz(-1.3335584) q[2];
sx q[2];
rz(-1.3875675) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(1.3685128) q[1];
sx q[1];
rz(-0.061731438) q[1];
sx q[1];
rz(1.3365937) q[1];
rz(-pi) q[2];
x q[2];
rz(-3.1114499) q[3];
sx q[3];
rz(-1.5557914) q[3];
sx q[3];
rz(-0.74149132) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(2.8162083) q[2];
sx q[2];
rz(-3.1367229) q[2];
sx q[2];
rz(-1.7083141) q[2];
rz(1.2895182) q[3];
sx q[3];
rz(-0.09434814) q[3];
sx q[3];
rz(-1.8452227) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0636487) q[0];
sx q[0];
rz(-2.1049121) q[0];
sx q[0];
rz(-2.5459469) q[0];
rz(-0.044364914) q[1];
sx q[1];
rz(-2.8652419) q[1];
sx q[1];
rz(-1.2038632) q[1];
rz(-0.2169256) q[2];
sx q[2];
rz(-1.1033162) q[2];
sx q[2];
rz(0.14324506) q[2];
rz(-1.4201319) q[3];
sx q[3];
rz(-0.41427783) q[3];
sx q[3];
rz(-2.8610341) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
