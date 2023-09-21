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
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.29205706) q[0];
sx q[0];
rz(-1.5229051) q[0];
sx q[0];
rz(3.1359262) q[0];
rz(-pi) q[1];
rz(-2.0179022) q[2];
sx q[2];
rz(-1.8748218) q[2];
sx q[2];
rz(-2.6927039) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(0.66558054) q[1];
sx q[1];
rz(-1.952938) q[1];
sx q[1];
rz(1.0784472) q[1];
rz(3.0311534) q[3];
sx q[3];
rz(-2.7613598) q[3];
sx q[3];
rz(-1.5766174) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(1.8346943) q[2];
sx q[2];
rz(-1.6480185) q[2];
sx q[2];
rz(-0.28960323) q[2];
rz(2.2662207) q[3];
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
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.61525476) q[0];
sx q[0];
rz(-2.2741788) q[0];
sx q[0];
rz(-0.34399024) q[0];
rz(-0.084331766) q[1];
sx q[1];
rz(-0.66939676) q[1];
sx q[1];
rz(1.7864236) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7536613) q[0];
sx q[0];
rz(-0.61457115) q[0];
sx q[0];
rz(-1.4580926) q[0];
rz(-pi) q[1];
x q[1];
rz(2.5002353) q[2];
sx q[2];
rz(-0.20198447) q[2];
sx q[2];
rz(-0.55822492) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(1.720571) q[1];
sx q[1];
rz(-0.45383006) q[1];
sx q[1];
rz(1.717091) q[1];
rz(-1.387272) q[3];
sx q[3];
rz(-1.4944544) q[3];
sx q[3];
rz(-0.21218382) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(2.8460059) q[2];
sx q[2];
rz(-0.83335525) q[2];
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
rz(pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.66353345) q[0];
sx q[0];
rz(-2.6038267) q[0];
sx q[0];
rz(-0.64087254) q[0];
rz(-2.3928941) q[1];
sx q[1];
rz(-2.1332032) q[1];
sx q[1];
rz(2.0764988) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.070806064) q[0];
sx q[0];
rz(-2.7822128) q[0];
sx q[0];
rz(1.552836) q[0];
rz(-pi) q[1];
rz(0.87419072) q[2];
sx q[2];
rz(-1.0993996) q[2];
sx q[2];
rz(3.0490321) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(1.1467495) q[1];
sx q[1];
rz(-0.61770505) q[1];
sx q[1];
rz(-1.5084933) q[1];
x q[2];
rz(-2.2933526) q[3];
sx q[3];
rz(-1.7602929) q[3];
sx q[3];
rz(-0.041785985) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(0.37725267) q[2];
sx q[2];
rz(-1.179402) q[2];
sx q[2];
rz(2.4242145) q[2];
rz(-0.68850368) q[3];
sx q[3];
rz(-2.513956) q[3];
sx q[3];
rz(2.8440516) q[3];
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
rz(-pi) q[0];
sx q[0];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3142969) q[0];
sx q[0];
rz(-1.2001487) q[0];
sx q[0];
rz(-2.8919019) q[0];
rz(-1.0149792) q[1];
sx q[1];
rz(-0.29622886) q[1];
sx q[1];
rz(0.011118523) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3664368) q[0];
sx q[0];
rz(-1.9754793) q[0];
sx q[0];
rz(2.0648271) q[0];
x q[1];
rz(0.87984933) q[2];
sx q[2];
rz(-0.86266154) q[2];
sx q[2];
rz(-2.326292) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-2.1383458) q[1];
sx q[1];
rz(-1.0838638) q[1];
sx q[1];
rz(-2.7702615) q[1];
x q[2];
rz(0.31828493) q[3];
sx q[3];
rz(-1.3878229) q[3];
sx q[3];
rz(1.1312248) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.49542385) q[2];
sx q[2];
rz(-1.255722) q[2];
sx q[2];
rz(-0.28309506) q[2];
rz(2.4781573) q[3];
sx q[3];
rz(-0.59316558) q[3];
sx q[3];
rz(-0.88808131) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0304612) q[0];
sx q[0];
rz(-0.24065329) q[0];
sx q[0];
rz(2.8097613) q[0];
rz(-2.6470673) q[1];
sx q[1];
rz(-1.8437513) q[1];
sx q[1];
rz(1.8146851) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.75816426) q[0];
sx q[0];
rz(-0.82793068) q[0];
sx q[0];
rz(-1.3616256) q[0];
rz(-0.97557108) q[2];
sx q[2];
rz(-2.274548) q[2];
sx q[2];
rz(-0.032035839) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-1.8718308) q[1];
sx q[1];
rz(-2.7360536) q[1];
sx q[1];
rz(-2.521442) q[1];
rz(-pi) q[2];
x q[2];
rz(1.9782449) q[3];
sx q[3];
rz(-1.3031928) q[3];
sx q[3];
rz(-1.0405376) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(1.1422687) q[2];
sx q[2];
rz(-0.48196718) q[2];
sx q[2];
rz(1.2456606) q[2];
rz(1.2549531) q[3];
sx q[3];
rz(-0.84147036) q[3];
sx q[3];
rz(0.83827034) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.44928837) q[0];
sx q[0];
rz(-0.0031539991) q[0];
sx q[0];
rz(-0.6814878) q[0];
rz(2.9340414) q[1];
sx q[1];
rz(-0.46336585) q[1];
sx q[1];
rz(-1.1157657) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.241248) q[0];
sx q[0];
rz(-1.1660518) q[0];
sx q[0];
rz(-0.83165283) q[0];
rz(-pi) q[1];
x q[1];
rz(2.7115466) q[2];
sx q[2];
rz(-0.94656813) q[2];
sx q[2];
rz(0.66832322) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(1.4286194) q[1];
sx q[1];
rz(-1.5884807) q[1];
sx q[1];
rz(2.5531205) q[1];
x q[2];
rz(2.4721018) q[3];
sx q[3];
rz(-1.4298425) q[3];
sx q[3];
rz(1.6657176) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-0.21268022) q[2];
sx q[2];
rz(-1.2968411) q[2];
sx q[2];
rz(0.35432717) q[2];
rz(2.8220693) q[3];
sx q[3];
rz(-1.1525681) q[3];
sx q[3];
rz(-2.7697146) q[3];
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
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0091244) q[0];
sx q[0];
rz(-2.9102944) q[0];
sx q[0];
rz(2.4672467) q[0];
rz(2.0293503) q[1];
sx q[1];
rz(-2.4770885) q[1];
sx q[1];
rz(2.5792714) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.092098504) q[0];
sx q[0];
rz(-1.5856165) q[0];
sx q[0];
rz(0.90256079) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.855905) q[2];
sx q[2];
rz(-0.34491587) q[2];
sx q[2];
rz(-3.0722741) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.14441227) q[1];
sx q[1];
rz(-1.4257396) q[1];
sx q[1];
rz(-0.1952862) q[1];
rz(-pi) q[2];
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
rz(1.101863) q[2];
sx q[2];
rz(-2.1477284) q[2];
sx q[2];
rz(-2.8015461) q[2];
rz(2.9240821) q[3];
sx q[3];
rz(-0.9483996) q[3];
sx q[3];
rz(-2.8229009) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6927004) q[0];
sx q[0];
rz(-3.0954439) q[0];
sx q[0];
rz(2.7451519) q[0];
rz(3.0026644) q[1];
sx q[1];
rz(-0.46008343) q[1];
sx q[1];
rz(-1.5213535) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8521261) q[0];
sx q[0];
rz(-2.9331101) q[0];
sx q[0];
rz(-2.8058488) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.0422158) q[2];
sx q[2];
rz(-1.2087012) q[2];
sx q[2];
rz(-0.21381703) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(0.24890451) q[1];
sx q[1];
rz(-2.3612594) q[1];
sx q[1];
rz(-0.08463879) q[1];
x q[2];
rz(0.89150724) q[3];
sx q[3];
rz(-1.4818947) q[3];
sx q[3];
rz(-0.31739435) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(0.56269318) q[2];
sx q[2];
rz(-1.0572628) q[2];
sx q[2];
rz(0.29433027) q[2];
rz(-2.0108022) q[3];
sx q[3];
rz(-1.7629938) q[3];
sx q[3];
rz(-2.14595) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[2];
sx q[3];
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
rz(2.6482553) q[0];
sx q[0];
rz(-2.2264037) q[0];
sx q[0];
rz(-0.35933581) q[0];
rz(-0.94611478) q[1];
sx q[1];
rz(-2.7455536) q[1];
sx q[1];
rz(-0.27063453) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.04979241) q[0];
sx q[0];
rz(-1.6066178) q[0];
sx q[0];
rz(-1.5396176) q[0];
rz(-pi) q[1];
rz(0.80910271) q[2];
sx q[2];
rz(-1.2247694) q[2];
sx q[2];
rz(-1.5835294) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.4295514) q[1];
sx q[1];
rz(-1.6152641) q[1];
sx q[1];
rz(3.0958423) q[1];
rz(-0.98468303) q[3];
sx q[3];
rz(-0.95041785) q[3];
sx q[3];
rz(-1.8970722) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(2.3140807) q[2];
sx q[2];
rz(-0.22958799) q[2];
sx q[2];
rz(0.45544004) q[2];
rz(2.3296302) q[3];
sx q[3];
rz(-0.88192314) q[3];
sx q[3];
rz(-2.5922095) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6311326) q[0];
sx q[0];
rz(-1.5083418) q[0];
sx q[0];
rz(0.73927885) q[0];
rz(-2.9108858) q[1];
sx q[1];
rz(-2.6735327) q[1];
sx q[1];
rz(-0.49490067) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4821645) q[0];
sx q[0];
rz(-0.48766252) q[0];
sx q[0];
rz(1.4898666) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.32785373) q[2];
sx q[2];
rz(-1.5640537) q[2];
sx q[2];
rz(-1.1525796) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(1.7086664) q[1];
sx q[1];
rz(-2.4927944) q[1];
sx q[1];
rz(1.385958) q[1];
rz(-1.6845735) q[3];
sx q[3];
rz(-2.4416231) q[3];
sx q[3];
rz(-1.1752807) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(0.13359244) q[2];
sx q[2];
rz(-1.091489) q[2];
sx q[2];
rz(2.8137394) q[2];
rz(-2.949529) q[3];
sx q[3];
rz(-2.8938507) q[3];
sx q[3];
rz(2.1081934) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1223758) q[0];
sx q[0];
rz(-1.5126956) q[0];
sx q[0];
rz(1.4737286) q[0];
rz(-2.3090251) q[1];
sx q[1];
rz(-1.7201798) q[1];
sx q[1];
rz(-1.7725772) q[1];
rz(1.0275502) q[2];
sx q[2];
rz(-1.3941358) q[2];
sx q[2];
rz(0.91640581) q[2];
rz(1.3151863) q[3];
sx q[3];
rz(-1.8288463) q[3];
sx q[3];
rz(2.0127206) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
