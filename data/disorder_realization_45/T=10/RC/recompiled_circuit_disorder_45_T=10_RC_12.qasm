OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-1.835445) q[0];
sx q[0];
rz(-0.68343502) q[0];
sx q[0];
rz(0.47877065) q[0];
rz(-3.1105644) q[1];
sx q[1];
rz(-1.9801158) q[1];
sx q[1];
rz(-0.64136139) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3135932) q[0];
sx q[0];
rz(-1.3324954) q[0];
sx q[0];
rz(0.39512623) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.1429206) q[2];
sx q[2];
rz(-1.2245721) q[2];
sx q[2];
rz(1.8510173) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(2.7176425) q[1];
sx q[1];
rz(-0.91361928) q[1];
sx q[1];
rz(2.1566725) q[1];
x q[2];
rz(1.1629421) q[3];
sx q[3];
rz(-2.1601094) q[3];
sx q[3];
rz(-2.4401963) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.59149867) q[2];
sx q[2];
rz(-1.2167565) q[2];
sx q[2];
rz(0.47810289) q[2];
rz(-1.452662) q[3];
sx q[3];
rz(-1.0457467) q[3];
sx q[3];
rz(2.5959192) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
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
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1801382) q[0];
sx q[0];
rz(-2.366876) q[0];
sx q[0];
rz(-2.1226728) q[0];
rz(1.6628751) q[1];
sx q[1];
rz(-2.5264085) q[1];
sx q[1];
rz(-2.5085124) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.86242005) q[0];
sx q[0];
rz(-2.3775568) q[0];
sx q[0];
rz(-3.0252302) q[0];
x q[1];
rz(0.74866809) q[2];
sx q[2];
rz(-2.0249172) q[2];
sx q[2];
rz(1.548896) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-0.1795579) q[1];
sx q[1];
rz(-1.9704559) q[1];
sx q[1];
rz(2.7707997) q[1];
rz(0.88926104) q[3];
sx q[3];
rz(-1.325843) q[3];
sx q[3];
rz(-0.8196866) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-2.7818266) q[2];
sx q[2];
rz(-0.40010139) q[2];
sx q[2];
rz(3.0241372) q[2];
rz(2.84058) q[3];
sx q[3];
rz(-1.8071226) q[3];
sx q[3];
rz(-1.3628179) q[3];
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
rz(-pi) q[0];
sx q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2837219) q[0];
sx q[0];
rz(-0.71325934) q[0];
sx q[0];
rz(0.088407956) q[0];
rz(-1.1075426) q[1];
sx q[1];
rz(-0.91845599) q[1];
sx q[1];
rz(2.450768) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4118977) q[0];
sx q[0];
rz(-1.6459961) q[0];
sx q[0];
rz(-0.15326432) q[0];
rz(-pi) q[1];
x q[1];
rz(2.8969953) q[2];
sx q[2];
rz(-2.3615712) q[2];
sx q[2];
rz(-0.56842677) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-1.1580116) q[1];
sx q[1];
rz(-2.9576655) q[1];
sx q[1];
rz(-1.4278825) q[1];
rz(-pi) q[2];
rz(0.37999837) q[3];
sx q[3];
rz(-1.899154) q[3];
sx q[3];
rz(-0.85193714) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-2.2174125) q[2];
sx q[2];
rz(-1.8841382) q[2];
sx q[2];
rz(3.0351191) q[2];
rz(1.7051833) q[3];
sx q[3];
rz(-2.7761716) q[3];
sx q[3];
rz(1.4829372) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0728264) q[0];
sx q[0];
rz(-2.9077353) q[0];
sx q[0];
rz(-2.6191214) q[0];
rz(-2.8126295) q[1];
sx q[1];
rz(-1.6429699) q[1];
sx q[1];
rz(0.55363384) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.75008167) q[0];
sx q[0];
rz(-0.9010074) q[0];
sx q[0];
rz(-2.1704587) q[0];
rz(-pi) q[1];
x q[1];
rz(0.81981084) q[2];
sx q[2];
rz(-2.1852583) q[2];
sx q[2];
rz(2.7053506) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-2.8071027) q[1];
sx q[1];
rz(-0.79542167) q[1];
sx q[1];
rz(2.9544178) q[1];
x q[2];
rz(-0.028285154) q[3];
sx q[3];
rz(-0.79704282) q[3];
sx q[3];
rz(-2.459211) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(2.5840977) q[2];
sx q[2];
rz(-2.2258874) q[2];
sx q[2];
rz(2.6468357) q[2];
rz(0.90302145) q[3];
sx q[3];
rz(-1.5938063) q[3];
sx q[3];
rz(-1.8516699) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6506127) q[0];
sx q[0];
rz(-0.74478331) q[0];
sx q[0];
rz(1.0850798) q[0];
rz(-2.1814573) q[1];
sx q[1];
rz(-1.1791869) q[1];
sx q[1];
rz(-0.18403149) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3071614) q[0];
sx q[0];
rz(-2.8330028) q[0];
sx q[0];
rz(-0.35085268) q[0];
rz(2.7765772) q[2];
sx q[2];
rz(-1.3114309) q[2];
sx q[2];
rz(0.58155453) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(1.4689456) q[1];
sx q[1];
rz(-0.99860672) q[1];
sx q[1];
rz(-2.1395626) q[1];
rz(2.9156978) q[3];
sx q[3];
rz(-2.4028006) q[3];
sx q[3];
rz(-1.6177288) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-2.2400143) q[2];
sx q[2];
rz(-2.7928536) q[2];
sx q[2];
rz(1.1408172) q[2];
rz(2.7187637) q[3];
sx q[3];
rz(-1.6641649) q[3];
sx q[3];
rz(2.0146577) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7891156) q[0];
sx q[0];
rz(-1.1081835) q[0];
sx q[0];
rz(-0.88678962) q[0];
rz(-2.9011762) q[1];
sx q[1];
rz(-1.0188894) q[1];
sx q[1];
rz(-2.9930847) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4473872) q[0];
sx q[0];
rz(-0.84852695) q[0];
sx q[0];
rz(1.3556051) q[0];
rz(-pi) q[1];
x q[1];
rz(0.44273419) q[2];
sx q[2];
rz(-0.25207439) q[2];
sx q[2];
rz(-0.67539757) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(1.8588184) q[1];
sx q[1];
rz(-2.8595279) q[1];
sx q[1];
rz(1.1689405) q[1];
x q[2];
rz(-1.968019) q[3];
sx q[3];
rz(-2.0052611) q[3];
sx q[3];
rz(1.43515) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(0.85577661) q[2];
sx q[2];
rz(-0.4824051) q[2];
sx q[2];
rz(-0.4883858) q[2];
rz(-2.6521902) q[3];
sx q[3];
rz(-0.92178744) q[3];
sx q[3];
rz(1.2667806) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4483036) q[0];
sx q[0];
rz(-1.332809) q[0];
sx q[0];
rz(0.81800246) q[0];
rz(1.8891634) q[1];
sx q[1];
rz(-2.7493582) q[1];
sx q[1];
rz(-2.0163527) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.97911191) q[0];
sx q[0];
rz(-1.9444124) q[0];
sx q[0];
rz(-3.0309249) q[0];
rz(-pi) q[1];
rz(2.0411885) q[2];
sx q[2];
rz(-2.8090968) q[2];
sx q[2];
rz(1.5770797) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-0.78912567) q[1];
sx q[1];
rz(-1.8450292) q[1];
sx q[1];
rz(1.2624361) q[1];
rz(-pi) q[2];
x q[2];
rz(2.4756487) q[3];
sx q[3];
rz(-2.0332554) q[3];
sx q[3];
rz(2.8921814) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-1.0730878) q[2];
sx q[2];
rz(-2.1606074) q[2];
sx q[2];
rz(0.64458624) q[2];
rz(-1.4792431) q[3];
sx q[3];
rz(-0.95696604) q[3];
sx q[3];
rz(-1.7361599) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.97061625) q[0];
sx q[0];
rz(-0.068844065) q[0];
sx q[0];
rz(-1.53565) q[0];
rz(-1.2212785) q[1];
sx q[1];
rz(-1.6284643) q[1];
sx q[1];
rz(2.1309526) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.37234026) q[0];
sx q[0];
rz(-1.4108037) q[0];
sx q[0];
rz(-0.82323797) q[0];
rz(-pi) q[1];
rz(2.0660731) q[2];
sx q[2];
rz(-1.3066548) q[2];
sx q[2];
rz(-1.8758945) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-2.752754) q[1];
sx q[1];
rz(-2.8201582) q[1];
sx q[1];
rz(2.3366117) q[1];
x q[2];
rz(-0.40057064) q[3];
sx q[3];
rz(-1.4860324) q[3];
sx q[3];
rz(0.29464196) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-1.5236139) q[2];
sx q[2];
rz(-2.8221059) q[2];
sx q[2];
rz(2.4475205) q[2];
rz(0.56898919) q[3];
sx q[3];
rz(-1.4469955) q[3];
sx q[3];
rz(0.93769658) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0555608) q[0];
sx q[0];
rz(-1.1431575) q[0];
sx q[0];
rz(-2.6468497) q[0];
rz(0.61839473) q[1];
sx q[1];
rz(-1.4952375) q[1];
sx q[1];
rz(-3.0659952) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8537124) q[0];
sx q[0];
rz(-1.259385) q[0];
sx q[0];
rz(-0.98543723) q[0];
rz(0.059968791) q[2];
sx q[2];
rz(-1.3955994) q[2];
sx q[2];
rz(2.2607468) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(1.1739993) q[1];
sx q[1];
rz(-0.6951957) q[1];
sx q[1];
rz(-2.980152) q[1];
x q[2];
rz(-1.3887614) q[3];
sx q[3];
rz(-0.958003) q[3];
sx q[3];
rz(0.36422563) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(0.33346924) q[2];
sx q[2];
rz(-1.7227017) q[2];
sx q[2];
rz(-1.3405651) q[2];
rz(0.30424413) q[3];
sx q[3];
rz(-2.2777568) q[3];
sx q[3];
rz(2.5568967) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
sx q[3];
rz(-pi) q[3];
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
rz(-1.8463523) q[0];
sx q[0];
rz(-1.3051935) q[0];
sx q[0];
rz(-2.9647968) q[0];
rz(1.8999752) q[1];
sx q[1];
rz(-2.5071564) q[1];
sx q[1];
rz(2.7005844) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.84687418) q[0];
sx q[0];
rz(-2.8537769) q[0];
sx q[0];
rz(-0.78886445) q[0];
x q[1];
rz(-1.4961692) q[2];
sx q[2];
rz(-1.7241038) q[2];
sx q[2];
rz(0.15664936) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.9778206) q[1];
sx q[1];
rz(-1.9107011) q[1];
sx q[1];
rz(2.0274859) q[1];
rz(-2.5173264) q[3];
sx q[3];
rz(-1.3253951) q[3];
sx q[3];
rz(2.0516968) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.5796154) q[2];
sx q[2];
rz(-2.5728971) q[2];
sx q[2];
rz(0.30187541) q[2];
rz(-0.89312303) q[3];
sx q[3];
rz(-1.8538657) q[3];
sx q[3];
rz(-0.20475234) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.72398913) q[0];
sx q[0];
rz(-1.2932734) q[0];
sx q[0];
rz(-1.4751157) q[0];
rz(-0.026731116) q[1];
sx q[1];
rz(-1.6865128) q[1];
sx q[1];
rz(-1.7105688) q[1];
rz(-2.2970207) q[2];
sx q[2];
rz(-1.9324586) q[2];
sx q[2];
rz(2.4697138) q[2];
rz(2.4206352) q[3];
sx q[3];
rz(-0.69098916) q[3];
sx q[3];
rz(0.50222764) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];