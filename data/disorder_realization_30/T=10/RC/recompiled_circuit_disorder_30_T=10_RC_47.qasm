OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(0.49139872) q[0];
sx q[0];
rz(2.8770652) q[0];
sx q[0];
rz(9.8192083) q[0];
rz(0.0061622942) q[1];
sx q[1];
rz(-0.34024629) q[1];
sx q[1];
rz(1.9415829) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8495356) q[0];
sx q[0];
rz(-1.6186876) q[0];
sx q[0];
rz(-0.0056664771) q[0];
rz(-pi) q[1];
rz(0.94304396) q[2];
sx q[2];
rz(-2.606751) q[2];
sx q[2];
rz(2.5778071) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-2.8438699) q[1];
sx q[1];
rz(-0.6134609) q[1];
sx q[1];
rz(-2.2754199) q[1];
rz(-pi) q[2];
rz(0.11043926) q[3];
sx q[3];
rz(-0.38023284) q[3];
sx q[3];
rz(-1.5766174) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(1.8346943) q[2];
sx q[2];
rz(-1.6480185) q[2];
sx q[2];
rz(-0.28960323) q[2];
rz(-0.87537193) q[3];
sx q[3];
rz(-2.1353728) q[3];
sx q[3];
rz(-3.0818821) q[3];
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
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.61525476) q[0];
sx q[0];
rz(-0.86741388) q[0];
sx q[0];
rz(-0.34399024) q[0];
rz(0.084331766) q[1];
sx q[1];
rz(-2.4721959) q[1];
sx q[1];
rz(1.7864236) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.25027572) q[0];
sx q[0];
rz(-2.180897) q[0];
sx q[0];
rz(-3.0623869) q[0];
rz(-pi) q[1];
x q[1];
rz(2.9789574) q[2];
sx q[2];
rz(-1.6911104) q[2];
sx q[2];
rz(-1.6441117) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-1.558074) q[1];
sx q[1];
rz(-1.1221702) q[1];
sx q[1];
rz(0.070986991) q[1];
rz(-pi) q[2];
rz(-1.387272) q[3];
sx q[3];
rz(-1.4944544) q[3];
sx q[3];
rz(2.9294088) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(0.29558674) q[2];
sx q[2];
rz(-2.3082374) q[2];
sx q[2];
rz(-2.6039092) q[2];
rz(2.6387571) q[3];
sx q[3];
rz(-0.42481315) q[3];
sx q[3];
rz(1.0104377) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.66353345) q[0];
sx q[0];
rz(-2.6038267) q[0];
sx q[0];
rz(-2.5007201) q[0];
rz(0.74869853) q[1];
sx q[1];
rz(-1.0083895) q[1];
sx q[1];
rz(-2.0764988) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.070806064) q[0];
sx q[0];
rz(-2.7822128) q[0];
sx q[0];
rz(-1.5887567) q[0];
x q[1];
rz(-2.5550585) q[2];
sx q[2];
rz(-0.9622935) q[2];
sx q[2];
rz(-1.3003132) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-2.0712191) q[1];
sx q[1];
rz(-2.1871236) q[1];
sx q[1];
rz(0.044205772) q[1];
rz(-pi) q[2];
rz(0.84824003) q[3];
sx q[3];
rz(-1.7602929) q[3];
sx q[3];
rz(3.0998067) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(2.76434) q[2];
sx q[2];
rz(-1.179402) q[2];
sx q[2];
rz(0.71737814) q[2];
rz(2.453089) q[3];
sx q[3];
rz(-2.513956) q[3];
sx q[3];
rz(2.8440516) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3142969) q[0];
sx q[0];
rz(-1.941444) q[0];
sx q[0];
rz(-2.8919019) q[0];
rz(-2.1266134) q[1];
sx q[1];
rz(-2.8453638) q[1];
sx q[1];
rz(-3.1304741) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3664368) q[0];
sx q[0];
rz(-1.1661134) q[0];
sx q[0];
rz(1.0767656) q[0];
rz(-pi) q[1];
rz(2.3035994) q[2];
sx q[2];
rz(-2.0760771) q[2];
sx q[2];
rz(2.879564) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(2.8324082) q[1];
sx q[1];
rz(-0.60317457) q[1];
sx q[1];
rz(-0.96997728) q[1];
rz(-1.3783781) q[3];
sx q[3];
rz(-1.8835861) q[3];
sx q[3];
rz(-0.49945143) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-0.49542385) q[2];
sx q[2];
rz(-1.8858706) q[2];
sx q[2];
rz(2.8584976) q[2];
rz(2.4781573) q[3];
sx q[3];
rz(-2.5484271) q[3];
sx q[3];
rz(-2.2535113) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
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
rz(3.0304612) q[0];
sx q[0];
rz(-2.9009394) q[0];
sx q[0];
rz(2.8097613) q[0];
rz(2.6470673) q[1];
sx q[1];
rz(-1.2978413) q[1];
sx q[1];
rz(1.8146851) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0622612) q[0];
sx q[0];
rz(-0.766303) q[0];
sx q[0];
rz(2.9192231) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.58381501) q[2];
sx q[2];
rz(-0.8875672) q[2];
sx q[2];
rz(-0.77606397) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-1.2111645) q[1];
sx q[1];
rz(-1.8976364) q[1];
sx q[1];
rz(-1.8153166) q[1];
rz(-0.96552403) q[3];
sx q[3];
rz(-11*pi/13) q[3];
sx q[3];
rz(-1.0799288) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.999324) q[2];
sx q[2];
rz(-2.6596255) q[2];
sx q[2];
rz(-1.8959321) q[2];
rz(1.8866395) q[3];
sx q[3];
rz(-2.3001223) q[3];
sx q[3];
rz(-2.3033223) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.44928837) q[0];
sx q[0];
rz(-0.0031539991) q[0];
sx q[0];
rz(-0.6814878) q[0];
rz(-0.20755126) q[1];
sx q[1];
rz(-0.46336585) q[1];
sx q[1];
rz(-1.1157657) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.67416699) q[0];
sx q[0];
rz(-0.90303991) q[0];
sx q[0];
rz(-0.52533124) q[0];
x q[1];
rz(0.43004604) q[2];
sx q[2];
rz(-0.94656813) q[2];
sx q[2];
rz(2.4732694) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-1.4286194) q[1];
sx q[1];
rz(-1.5884807) q[1];
sx q[1];
rz(-2.5531205) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.7498155) q[3];
sx q[3];
rz(-2.2324623) q[3];
sx q[3];
rz(0.20565198) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.9289124) q[2];
sx q[2];
rz(-1.2968411) q[2];
sx q[2];
rz(2.7872655) q[2];
rz(2.8220693) q[3];
sx q[3];
rz(-1.1525681) q[3];
sx q[3];
rz(-2.7697146) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
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
rz(-0.56232125) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0494941) q[0];
sx q[0];
rz(-1.5559762) q[0];
sx q[0];
rz(-2.2390319) q[0];
rz(2.8096335) q[2];
sx q[2];
rz(-1.6662285) q[2];
sx q[2];
rz(-1.7711668) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-2.9971804) q[1];
sx q[1];
rz(-1.715853) q[1];
sx q[1];
rz(-0.1952862) q[1];
x q[2];
rz(-1.4427156) q[3];
sx q[3];
rz(-2.5411385) q[3];
sx q[3];
rz(-1.2680935) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.0397296) q[2];
sx q[2];
rz(-0.9938643) q[2];
sx q[2];
rz(2.8015461) q[2];
rz(0.21751054) q[3];
sx q[3];
rz(-0.9483996) q[3];
sx q[3];
rz(-0.31869179) q[3];
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
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6927004) q[0];
sx q[0];
rz(-0.046148766) q[0];
sx q[0];
rz(2.7451519) q[0];
rz(-0.13892826) q[1];
sx q[1];
rz(-2.6815092) q[1];
sx q[1];
rz(-1.6202392) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0884468) q[0];
sx q[0];
rz(-1.3741115) q[0];
sx q[0];
rz(-1.6403857) q[0];
rz(-0.87562008) q[2];
sx q[2];
rz(-2.5555829) q[2];
sx q[2];
rz(-2.3919174) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-1.382114) q[1];
sx q[1];
rz(-1.5112875) q[1];
sx q[1];
rz(-2.3630523) q[1];
rz(-pi) q[2];
rz(-0.89150724) q[3];
sx q[3];
rz(-1.4818947) q[3];
sx q[3];
rz(0.31739435) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.5788995) q[2];
sx q[2];
rz(-1.0572628) q[2];
sx q[2];
rz(-2.8472624) q[2];
rz(2.0108022) q[3];
sx q[3];
rz(-1.3785988) q[3];
sx q[3];
rz(-2.14595) q[3];
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
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6482553) q[0];
sx q[0];
rz(-0.91518891) q[0];
sx q[0];
rz(-0.35933581) q[0];
rz(2.1954779) q[1];
sx q[1];
rz(-0.39603907) q[1];
sx q[1];
rz(0.27063453) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.7662738) q[0];
sx q[0];
rz(-0.047485504) q[0];
sx q[0];
rz(-0.71592285) q[0];
rz(-pi) q[1];
x q[1];
rz(1.0893906) q[2];
sx q[2];
rz(-0.82197661) q[2];
sx q[2];
rz(2.7871728) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-2.4295514) q[1];
sx q[1];
rz(-1.5263285) q[1];
sx q[1];
rz(0.045750381) q[1];
rz(2.4828033) q[3];
sx q[3];
rz(-0.82595347) q[3];
sx q[3];
rz(0.39289075) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(2.3140807) q[2];
sx q[2];
rz(-2.9120047) q[2];
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
rz(pi/2) q[1];
rz(-pi) q[2];
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
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.51046002) q[0];
sx q[0];
rz(-1.5083418) q[0];
sx q[0];
rz(-0.73927885) q[0];
rz(0.23070681) q[1];
sx q[1];
rz(-0.46805996) q[1];
sx q[1];
rz(-2.646692) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.98289821) q[0];
sx q[0];
rz(-1.532908) q[0];
sx q[0];
rz(1.0844896) q[0];
rz(-pi) q[1];
rz(-0.020936326) q[2];
sx q[2];
rz(-2.8136721) q[2];
sx q[2];
rz(-0.39839572) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-1.478185) q[1];
sx q[1];
rz(-2.2067398) q[1];
sx q[1];
rz(-0.138476) q[1];
rz(-pi) q[2];
rz(-0.095330843) q[3];
sx q[3];
rz(-2.2653326) q[3];
sx q[3];
rz(1.3235843) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-3.0080002) q[2];
sx q[2];
rz(-2.0501037) q[2];
sx q[2];
rz(-2.8137394) q[2];
rz(-0.19206364) q[3];
sx q[3];
rz(-2.8938507) q[3];
sx q[3];
rz(1.0333992) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1223758) q[0];
sx q[0];
rz(-1.5126956) q[0];
sx q[0];
rz(1.4737286) q[0];
rz(-0.83256759) q[1];
sx q[1];
rz(-1.4214129) q[1];
sx q[1];
rz(1.3690154) q[1];
rz(-2.1140425) q[2];
sx q[2];
rz(-1.3941358) q[2];
sx q[2];
rz(0.91640581) q[2];
rz(1.8264063) q[3];
sx q[3];
rz(-1.3127463) q[3];
sx q[3];
rz(-1.128872) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
