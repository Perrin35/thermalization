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
rz(-1.5440829) q[0];
sx q[0];
rz(-1.7680327) q[0];
sx q[0];
rz(-1.0184259) q[0];
rz(-0.51813689) q[1];
sx q[1];
rz(-2.3845446) q[1];
sx q[1];
rz(-2.5098324) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3903094) q[0];
sx q[0];
rz(-1.3645118) q[0];
sx q[0];
rz(0.12796107) q[0];
rz(-pi) q[1];
rz(-0.024876923) q[2];
sx q[2];
rz(-1.2681677) q[2];
sx q[2];
rz(-0.51371117) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(2.2871197) q[1];
sx q[1];
rz(-1.8504256) q[1];
sx q[1];
rz(0.92052144) q[1];
rz(-1.8970892) q[3];
sx q[3];
rz(-1.4832116) q[3];
sx q[3];
rz(2.0503941) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.6112001) q[2];
sx q[2];
rz(-2.7860614) q[2];
sx q[2];
rz(-1.6543039) q[2];
rz(-2.2330331) q[3];
sx q[3];
rz(-0.23659758) q[3];
sx q[3];
rz(-1.0049741) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0205883) q[0];
sx q[0];
rz(-1.7089184) q[0];
sx q[0];
rz(-0.16227907) q[0];
rz(-0.24457112) q[1];
sx q[1];
rz(-1.9385612) q[1];
sx q[1];
rz(2.8499106) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8732018) q[0];
sx q[0];
rz(-1.8531728) q[0];
sx q[0];
rz(1.3703521) q[0];
rz(-pi) q[1];
x q[1];
rz(1.2873093) q[2];
sx q[2];
rz(-1.3559196) q[2];
sx q[2];
rz(2.7708997) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-1.6704935) q[1];
sx q[1];
rz(-2.3355977) q[1];
sx q[1];
rz(2.8189895) q[1];
x q[2];
rz(-0.41145153) q[3];
sx q[3];
rz(-1.1845767) q[3];
sx q[3];
rz(-2.7038108) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-0.67819277) q[2];
sx q[2];
rz(-2.2726111) q[2];
sx q[2];
rz(2.0657067) q[2];
rz(1.894527) q[3];
sx q[3];
rz(-1.5445292) q[3];
sx q[3];
rz(1.6943078) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4397044) q[0];
sx q[0];
rz(-2.0151558) q[0];
sx q[0];
rz(2.9578748) q[0];
rz(-1.6366929) q[1];
sx q[1];
rz(-0.74438649) q[1];
sx q[1];
rz(-2.9768129) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4713584) q[0];
sx q[0];
rz(-2.5557098) q[0];
sx q[0];
rz(1.5536932) q[0];
rz(-0.38305958) q[2];
sx q[2];
rz(-2.428741) q[2];
sx q[2];
rz(2.0682316) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(0.91080571) q[1];
sx q[1];
rz(-1.1890829) q[1];
sx q[1];
rz(2.4934105) q[1];
rz(-pi) q[2];
rz(-2.1117049) q[3];
sx q[3];
rz(-0.60324429) q[3];
sx q[3];
rz(1.3689976) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-3.034721) q[2];
sx q[2];
rz(-1.9811337) q[2];
sx q[2];
rz(-1.1064233) q[2];
rz(2.4404081) q[3];
sx q[3];
rz(-1.8132352) q[3];
sx q[3];
rz(0.4655233) q[3];
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
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.768854) q[0];
sx q[0];
rz(-1.1356069) q[0];
sx q[0];
rz(2.1970774) q[0];
rz(1.6230029) q[1];
sx q[1];
rz(-0.98097643) q[1];
sx q[1];
rz(1.6015582) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.9409143) q[0];
sx q[0];
rz(-0.8884065) q[0];
sx q[0];
rz(2.3290578) q[0];
rz(-pi) q[1];
rz(0.77727274) q[2];
sx q[2];
rz(-0.16675719) q[2];
sx q[2];
rz(0.47698944) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-0.67805144) q[1];
sx q[1];
rz(-1.9423331) q[1];
sx q[1];
rz(2.3080565) q[1];
x q[2];
rz(0.047961162) q[3];
sx q[3];
rz(-0.75100079) q[3];
sx q[3];
rz(2.6329071) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(3.1314142) q[2];
sx q[2];
rz(-0.62128908) q[2];
sx q[2];
rz(0.16769257) q[2];
rz(-0.0532648) q[3];
sx q[3];
rz(-2.0361418) q[3];
sx q[3];
rz(1.3767327) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5662956) q[0];
sx q[0];
rz(-3.0360041) q[0];
sx q[0];
rz(-2.8422624) q[0];
rz(0.37295595) q[1];
sx q[1];
rz(-1.2495709) q[1];
sx q[1];
rz(1.4704871) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2946314) q[0];
sx q[0];
rz(-1.9982013) q[0];
sx q[0];
rz(0.60890108) q[0];
x q[1];
rz(-2.9048033) q[2];
sx q[2];
rz(-2.6663187) q[2];
sx q[2];
rz(-1.4314112) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.59488397) q[1];
sx q[1];
rz(-2.1404033) q[1];
sx q[1];
rz(-0.96086603) q[1];
x q[2];
rz(0.94909747) q[3];
sx q[3];
rz(-1.4273564) q[3];
sx q[3];
rz(2.6196817) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-1.9386998) q[2];
sx q[2];
rz(-2.4788269) q[2];
sx q[2];
rz(-2.0311484) q[2];
rz(0.86197305) q[3];
sx q[3];
rz(-1.5098666) q[3];
sx q[3];
rz(1.8814258) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.92457572) q[0];
sx q[0];
rz(-0.68704263) q[0];
sx q[0];
rz(0.4050912) q[0];
rz(0.336126) q[1];
sx q[1];
rz(-1.4210217) q[1];
sx q[1];
rz(0.95132557) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8444871) q[0];
sx q[0];
rz(-1.9713624) q[0];
sx q[0];
rz(1.8482918) q[0];
x q[1];
rz(0.73619618) q[2];
sx q[2];
rz(-1.7288107) q[2];
sx q[2];
rz(0.73537961) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(1.8321633) q[1];
sx q[1];
rz(-2.1346666) q[1];
sx q[1];
rz(-1.840074) q[1];
rz(-0.59744617) q[3];
sx q[3];
rz(-2.8077112) q[3];
sx q[3];
rz(2.2528439) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-2.090791) q[2];
sx q[2];
rz(-1.7369221) q[2];
sx q[2];
rz(2.9108099) q[2];
rz(0.18203059) q[3];
sx q[3];
rz(-2.5735276) q[3];
sx q[3];
rz(0.52869421) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
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
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.35293216) q[0];
sx q[0];
rz(-3.1020628) q[0];
sx q[0];
rz(-0.93210644) q[0];
rz(2.508029) q[1];
sx q[1];
rz(-2.0220058) q[1];
sx q[1];
rz(1.8036141) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.47565002) q[0];
sx q[0];
rz(-1.4650657) q[0];
sx q[0];
rz(-0.08203489) q[0];
rz(-0.39057486) q[2];
sx q[2];
rz(-2.8468067) q[2];
sx q[2];
rz(-0.24402555) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(1.6465969) q[1];
sx q[1];
rz(-0.75237583) q[1];
sx q[1];
rz(-0.97480358) q[1];
rz(-pi) q[2];
rz(-0.49832817) q[3];
sx q[3];
rz(-1.5987724) q[3];
sx q[3];
rz(1.2714112) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-0.46447095) q[2];
sx q[2];
rz(-1.9147583) q[2];
sx q[2];
rz(-1.202549) q[2];
rz(2.7770212) q[3];
sx q[3];
rz(-2.4489844) q[3];
sx q[3];
rz(-2.306849) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6637591) q[0];
sx q[0];
rz(-2.5680225) q[0];
sx q[0];
rz(-2.9649576) q[0];
rz(-0.44003507) q[1];
sx q[1];
rz(-1.1853848) q[1];
sx q[1];
rz(-2.1766591) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2595183) q[0];
sx q[0];
rz(-1.4384369) q[0];
sx q[0];
rz(-2.9455723) q[0];
rz(-pi) q[1];
x q[1];
rz(0.60243269) q[2];
sx q[2];
rz(-2.0332073) q[2];
sx q[2];
rz(1.3585857) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(1.392004) q[1];
sx q[1];
rz(-0.67396213) q[1];
sx q[1];
rz(-2.1753009) q[1];
rz(-pi) q[2];
rz(1.8252402) q[3];
sx q[3];
rz(-1.4973876) q[3];
sx q[3];
rz(2.0919098) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-2.9289916) q[2];
sx q[2];
rz(-1.5233728) q[2];
sx q[2];
rz(3.1316481) q[2];
rz(3.0250004) q[3];
sx q[3];
rz(-0.3392342) q[3];
sx q[3];
rz(2.5998083) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.56124878) q[0];
sx q[0];
rz(-1.2224226) q[0];
sx q[0];
rz(-2.5196581) q[0];
rz(1.7033345) q[1];
sx q[1];
rz(-0.57241264) q[1];
sx q[1];
rz(0.66351801) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.91231649) q[0];
sx q[0];
rz(-1.7890507) q[0];
sx q[0];
rz(2.403232) q[0];
x q[1];
rz(-1.8542669) q[2];
sx q[2];
rz(-1.4525692) q[2];
sx q[2];
rz(-0.17519874) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-0.27706669) q[1];
sx q[1];
rz(-2.4063595) q[1];
sx q[1];
rz(1.0864054) q[1];
rz(2.8264753) q[3];
sx q[3];
rz(-2.8459918) q[3];
sx q[3];
rz(2.78418) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-0.62349391) q[2];
sx q[2];
rz(-1.7524717) q[2];
sx q[2];
rz(2.7790879) q[2];
rz(-2.6643961) q[3];
sx q[3];
rz(-1.0537182) q[3];
sx q[3];
rz(1.1227192) q[3];
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
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0839194) q[0];
sx q[0];
rz(-2.4102983) q[0];
sx q[0];
rz(-2.2398563) q[0];
rz(-2.4665191) q[1];
sx q[1];
rz(-2.156064) q[1];
sx q[1];
rz(2.9973082) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0269723) q[0];
sx q[0];
rz(-1.2204224) q[0];
sx q[0];
rz(3.0536985) q[0];
rz(0.21964964) q[2];
sx q[2];
rz(-1.4712508) q[2];
sx q[2];
rz(1.8451898) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(0.044202494) q[1];
sx q[1];
rz(-2.4501738) q[1];
sx q[1];
rz(-2.7425062) q[1];
rz(-1.3809526) q[3];
sx q[3];
rz(-1.1355054) q[3];
sx q[3];
rz(1.3224885) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.5054063) q[2];
sx q[2];
rz(-1.9289086) q[2];
sx q[2];
rz(2.9409161) q[2];
rz(2.0319273) q[3];
sx q[3];
rz(-0.4626914) q[3];
sx q[3];
rz(2.5698575) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
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
rz(1.5526445) q[0];
sx q[0];
rz(-0.968796) q[0];
sx q[0];
rz(2.0198685) q[0];
rz(0.48925346) q[1];
sx q[1];
rz(-1.4514634) q[1];
sx q[1];
rz(-1.0101752) q[1];
rz(1.3623753) q[2];
sx q[2];
rz(-0.42283146) q[2];
sx q[2];
rz(-2.6953807) q[2];
rz(1.5420001) q[3];
sx q[3];
rz(-1.1287108) q[3];
sx q[3];
rz(-2.5870833) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
