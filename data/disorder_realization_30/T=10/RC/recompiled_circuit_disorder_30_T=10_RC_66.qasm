OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-2.6501939) q[0];
sx q[0];
rz(-2.8770652) q[0];
sx q[0];
rz(-2.7471623) q[0];
rz(0.0061622942) q[1];
sx q[1];
rz(2.8013464) q[1];
sx q[1];
rz(10.624788) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.29205706) q[0];
sx q[0];
rz(-1.5229051) q[0];
sx q[0];
rz(-0.0056664771) q[0];
rz(-0.33485246) q[2];
sx q[2];
rz(-1.1455673) q[2];
sx q[2];
rz(1.877117) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(0.66558054) q[1];
sx q[1];
rz(-1.1886547) q[1];
sx q[1];
rz(-1.0784472) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.5267738) q[3];
sx q[3];
rz(-1.9485954) q[3];
sx q[3];
rz(1.4461185) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(1.3068984) q[2];
sx q[2];
rz(-1.4935741) q[2];
sx q[2];
rz(2.8519894) q[2];
rz(-0.87537193) q[3];
sx q[3];
rz(-2.1353728) q[3];
sx q[3];
rz(0.059710596) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.61525476) q[0];
sx q[0];
rz(-0.86741388) q[0];
sx q[0];
rz(0.34399024) q[0];
rz(0.084331766) q[1];
sx q[1];
rz(-0.66939676) q[1];
sx q[1];
rz(1.3551691) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8913169) q[0];
sx q[0];
rz(-2.180897) q[0];
sx q[0];
rz(-3.0623869) q[0];
x q[1];
rz(-1.4488892) q[2];
sx q[2];
rz(-1.4093471) q[2];
sx q[2];
rz(0.09300692) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-1.720571) q[1];
sx q[1];
rz(-0.45383006) q[1];
sx q[1];
rz(-1.717091) q[1];
rz(-3.063952) q[3];
sx q[3];
rz(-1.3878126) q[3];
sx q[3];
rz(-1.3444572) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(0.29558674) q[2];
sx q[2];
rz(-0.83335525) q[2];
sx q[2];
rz(2.6039092) q[2];
rz(0.50283557) q[3];
sx q[3];
rz(-0.42481315) q[3];
sx q[3];
rz(-1.0104377) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.66353345) q[0];
sx q[0];
rz(-2.6038267) q[0];
sx q[0];
rz(-2.5007201) q[0];
rz(2.3928941) q[1];
sx q[1];
rz(-2.1332032) q[1];
sx q[1];
rz(1.0650939) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0707866) q[0];
sx q[0];
rz(-0.35937989) q[0];
sx q[0];
rz(1.5887567) q[0];
rz(-pi) q[1];
rz(-0.87419072) q[2];
sx q[2];
rz(-1.0993996) q[2];
sx q[2];
rz(0.092560571) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-2.0712191) q[1];
sx q[1];
rz(-0.95446903) q[1];
sx q[1];
rz(-0.044205772) q[1];
x q[2];
rz(-2.2933526) q[3];
sx q[3];
rz(-1.3812997) q[3];
sx q[3];
rz(-3.0998067) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-2.76434) q[2];
sx q[2];
rz(-1.9621907) q[2];
sx q[2];
rz(-2.4242145) q[2];
rz(-2.453089) q[3];
sx q[3];
rz(-0.62763667) q[3];
sx q[3];
rz(2.8440516) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
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
rz(2.3142969) q[0];
sx q[0];
rz(-1.941444) q[0];
sx q[0];
rz(-0.24969077) q[0];
rz(2.1266134) q[1];
sx q[1];
rz(-0.29622886) q[1];
sx q[1];
rz(0.011118523) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.5866833) q[0];
sx q[0];
rz(-2.0218098) q[0];
sx q[0];
rz(0.45278544) q[0];
x q[1];
rz(2.2617433) q[2];
sx q[2];
rz(-0.86266154) q[2];
sx q[2];
rz(2.326292) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(2.8324082) q[1];
sx q[1];
rz(-2.5384181) q[1];
sx q[1];
rz(-2.1716154) q[1];
rz(-0.31828493) q[3];
sx q[3];
rz(-1.3878229) q[3];
sx q[3];
rz(2.0103679) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(2.6461688) q[2];
sx q[2];
rz(-1.255722) q[2];
sx q[2];
rz(0.28309506) q[2];
rz(0.66343534) q[3];
sx q[3];
rz(-0.59316558) q[3];
sx q[3];
rz(-2.2535113) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.11113142) q[0];
sx q[0];
rz(-2.9009394) q[0];
sx q[0];
rz(2.8097613) q[0];
rz(0.49452531) q[1];
sx q[1];
rz(-1.2978413) q[1];
sx q[1];
rz(-1.8146851) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.6700232) q[0];
sx q[0];
rz(-1.7243392) q[0];
sx q[0];
rz(-0.75385401) q[0];
x q[1];
rz(2.5577776) q[2];
sx q[2];
rz(-0.8875672) q[2];
sx q[2];
rz(-0.77606397) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(1.8718308) q[1];
sx q[1];
rz(-2.7360536) q[1];
sx q[1];
rz(2.521442) q[1];
x q[2];
rz(-2.1760686) q[3];
sx q[3];
rz(-11*pi/13) q[3];
sx q[3];
rz(-2.0616639) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(1.999324) q[2];
sx q[2];
rz(-2.6596255) q[2];
sx q[2];
rz(1.2456606) q[2];
rz(-1.8866395) q[3];
sx q[3];
rz(-0.84147036) q[3];
sx q[3];
rz(-2.3033223) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6923043) q[0];
sx q[0];
rz(-0.0031539991) q[0];
sx q[0];
rz(-0.6814878) q[0];
rz(0.20755126) q[1];
sx q[1];
rz(-2.6782268) q[1];
sx q[1];
rz(2.025827) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.078243144) q[0];
sx q[0];
rz(-2.3176498) q[0];
sx q[0];
rz(-2.1372165) q[0];
rz(-pi) q[1];
rz(-2.240928) q[2];
sx q[2];
rz(-1.2256983) q[2];
sx q[2];
rz(-2.5010441) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(1.4286194) q[1];
sx q[1];
rz(-1.5884807) q[1];
sx q[1];
rz(-0.58847217) q[1];
rz(-pi) q[2];
x q[2];
rz(1.3917771) q[3];
sx q[3];
rz(-2.2324623) q[3];
sx q[3];
rz(0.20565198) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-2.9289124) q[2];
sx q[2];
rz(-1.2968411) q[2];
sx q[2];
rz(2.7872655) q[2];
rz(-0.3195233) q[3];
sx q[3];
rz(-1.9890246) q[3];
sx q[3];
rz(2.7697146) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
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
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1324683) q[0];
sx q[0];
rz(-2.9102944) q[0];
sx q[0];
rz(-2.4672467) q[0];
rz(2.0293503) q[1];
sx q[1];
rz(-2.4770885) q[1];
sx q[1];
rz(2.5792714) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0494941) q[0];
sx q[0];
rz(-1.5559762) q[0];
sx q[0];
rz(0.90256079) q[0];
x q[1];
rz(0.33195915) q[2];
sx q[2];
rz(-1.6662285) q[2];
sx q[2];
rz(-1.3704259) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(0.14441227) q[1];
sx q[1];
rz(-1.715853) q[1];
sx q[1];
rz(2.9463065) q[1];
rz(3.0543442) q[3];
sx q[3];
rz(-0.97594075) q[3];
sx q[3];
rz(-2.0283386) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-1.101863) q[2];
sx q[2];
rz(-2.1477284) q[2];
sx q[2];
rz(-0.34004655) q[2];
rz(-0.21751054) q[3];
sx q[3];
rz(-0.9483996) q[3];
sx q[3];
rz(-2.8229009) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
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
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.44889221) q[0];
sx q[0];
rz(-0.046148766) q[0];
sx q[0];
rz(2.7451519) q[0];
rz(-0.13892826) q[1];
sx q[1];
rz(-2.6815092) q[1];
sx q[1];
rz(1.5213535) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6103219) q[0];
sx q[0];
rz(-1.5025508) q[0];
sx q[0];
rz(-2.9444429) q[0];
rz(-0.87562008) q[2];
sx q[2];
rz(-2.5555829) q[2];
sx q[2];
rz(-2.3919174) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-3.0115067) q[1];
sx q[1];
rz(-2.3475921) q[1];
sx q[1];
rz(1.6542875) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.2500854) q[3];
sx q[3];
rz(-1.659698) q[3];
sx q[3];
rz(0.31739435) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(2.5788995) q[2];
sx q[2];
rz(-1.0572628) q[2];
sx q[2];
rz(2.8472624) q[2];
rz(-1.1307905) q[3];
sx q[3];
rz(-1.7629938) q[3];
sx q[3];
rz(-0.99564266) q[3];
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
x q[0];
rz(-pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.49333736) q[0];
sx q[0];
rz(-2.2264037) q[0];
sx q[0];
rz(2.7822568) q[0];
rz(-0.94611478) q[1];
sx q[1];
rz(-0.39603907) q[1];
sx q[1];
rz(0.27063453) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.7662738) q[0];
sx q[0];
rz(-0.047485504) q[0];
sx q[0];
rz(2.4256698) q[0];
x q[1];
rz(2.3324899) q[2];
sx q[2];
rz(-1.9168233) q[2];
sx q[2];
rz(1.5580633) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(3.0535069) q[1];
sx q[1];
rz(-0.063789531) q[1];
sx q[1];
rz(0.77168685) q[1];
x q[2];
rz(2.1569096) q[3];
sx q[3];
rz(-0.95041785) q[3];
sx q[3];
rz(1.2445205) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(0.82751194) q[2];
sx q[2];
rz(-0.22958799) q[2];
sx q[2];
rz(-0.45544004) q[2];
rz(2.3296302) q[3];
sx q[3];
rz(-2.2596695) q[3];
sx q[3];
rz(-0.5493831) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6311326) q[0];
sx q[0];
rz(-1.5083418) q[0];
sx q[0];
rz(2.4023138) q[0];
rz(-0.23070681) q[1];
sx q[1];
rz(-0.46805996) q[1];
sx q[1];
rz(2.646692) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5737168) q[0];
sx q[0];
rz(-1.084869) q[0];
sx q[0];
rz(3.0987415) q[0];
x q[1];
rz(2.8137389) q[2];
sx q[2];
rz(-1.5640537) q[2];
sx q[2];
rz(-1.1525796) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.6634076) q[1];
sx q[1];
rz(-0.93485281) q[1];
sx q[1];
rz(3.0031167) q[1];
x q[2];
rz(-1.4570191) q[3];
sx q[3];
rz(-0.69996951) q[3];
sx q[3];
rz(-1.1752807) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(3.0080002) q[2];
sx q[2];
rz(-2.0501037) q[2];
sx q[2];
rz(-0.32785329) q[2];
rz(0.19206364) q[3];
sx q[3];
rz(-0.24774194) q[3];
sx q[3];
rz(1.0333992) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
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
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1223758) q[0];
sx q[0];
rz(-1.6288971) q[0];
sx q[0];
rz(-1.6678641) q[0];
rz(-2.3090251) q[1];
sx q[1];
rz(-1.7201798) q[1];
sx q[1];
rz(-1.7725772) q[1];
rz(-1.9033296) q[2];
sx q[2];
rz(-0.56849545) q[2];
sx q[2];
rz(-0.93760437) q[2];
rz(-2.3776688) q[3];
sx q[3];
rz(-2.7803956) q[3];
sx q[3];
rz(-1.9261388) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];