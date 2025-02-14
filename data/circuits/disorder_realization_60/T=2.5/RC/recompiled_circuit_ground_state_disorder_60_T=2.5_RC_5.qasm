OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(0.085091703) q[0];
sx q[0];
rz(-0.29274517) q[0];
sx q[0];
rz(-0.91896397) q[0];
rz(-0.38209823) q[1];
sx q[1];
rz(-2.9736019) q[1];
sx q[1];
rz(1.1408495) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9716233) q[0];
sx q[0];
rz(-0.21339082) q[0];
sx q[0];
rz(2.1687228) q[0];
rz(-pi) q[1];
rz(1.2873285) q[2];
sx q[2];
rz(-1.8745443) q[2];
sx q[2];
rz(1.8423353) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(0.39841043) q[1];
sx q[1];
rz(-2.7355621) q[1];
sx q[1];
rz(-0.34161411) q[1];
rz(-pi) q[2];
rz(-1.857418) q[3];
sx q[3];
rz(-2.686073) q[3];
sx q[3];
rz(0.16715967) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(0.69316489) q[2];
sx q[2];
rz(-1.0591256) q[2];
sx q[2];
rz(-1.1761752) q[2];
rz(3.0991992) q[3];
sx q[3];
rz(-1.8040801) q[3];
sx q[3];
rz(0.5347518) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6723044) q[0];
sx q[0];
rz(-2.7870218) q[0];
sx q[0];
rz(2.9571423) q[0];
rz(-0.09672673) q[1];
sx q[1];
rz(-1.6177142) q[1];
sx q[1];
rz(1.2299889) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7589129) q[0];
sx q[0];
rz(-1.2059187) q[0];
sx q[0];
rz(-1.7880102) q[0];
x q[1];
rz(-1.5225836) q[2];
sx q[2];
rz(-2.1297751) q[2];
sx q[2];
rz(-0.79254442) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(1.1538196) q[1];
sx q[1];
rz(-1.7708988) q[1];
sx q[1];
rz(0.2421744) q[1];
rz(-pi) q[2];
rz(1.6685247) q[3];
sx q[3];
rz(-2.7616581) q[3];
sx q[3];
rz(-1.987628) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-2.7868598) q[2];
sx q[2];
rz(-0.41149461) q[2];
sx q[2];
rz(0.011431781) q[2];
rz(1.8276021) q[3];
sx q[3];
rz(-1.3649789) q[3];
sx q[3];
rz(-2.438681) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
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
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3979724) q[0];
sx q[0];
rz(-3.008606) q[0];
sx q[0];
rz(-2.5308894) q[0];
rz(0.01344219) q[1];
sx q[1];
rz(-1.8553269) q[1];
sx q[1];
rz(0.59649831) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0666351) q[0];
sx q[0];
rz(-2.5078815) q[0];
sx q[0];
rz(1.8541965) q[0];
x q[1];
rz(0.99278583) q[2];
sx q[2];
rz(-1.3137378) q[2];
sx q[2];
rz(1.3959194) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(2.8104659) q[1];
sx q[1];
rz(-0.8672204) q[1];
sx q[1];
rz(-1.1570279) q[1];
rz(1.0168996) q[3];
sx q[3];
rz(-1.5988598) q[3];
sx q[3];
rz(1.156776) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-0.67768031) q[2];
sx q[2];
rz(-1.3119421) q[2];
sx q[2];
rz(-0.00027351969) q[2];
rz(1.6655946) q[3];
sx q[3];
rz(-0.6690343) q[3];
sx q[3];
rz(-3.0345548) q[3];
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
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.45604712) q[0];
sx q[0];
rz(-0.16655971) q[0];
sx q[0];
rz(-1.9786932) q[0];
rz(-2.1641425) q[1];
sx q[1];
rz(-1.6012499) q[1];
sx q[1];
rz(-1.5862484) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.68846164) q[0];
sx q[0];
rz(-1.5630088) q[0];
sx q[0];
rz(-1.5611737) q[0];
rz(3.0968148) q[2];
sx q[2];
rz(-2.3557202) q[2];
sx q[2];
rz(-1.7742243) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(1.697243) q[1];
sx q[1];
rz(-1.0740136) q[1];
sx q[1];
rz(-1.5543943) q[1];
x q[2];
rz(-1.076872) q[3];
sx q[3];
rz(-1.3528429) q[3];
sx q[3];
rz(2.2214784) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-2.1996475) q[2];
sx q[2];
rz(-2.6298099) q[2];
sx q[2];
rz(0.23435782) q[2];
rz(-0.50218454) q[3];
sx q[3];
rz(-2.1000704) q[3];
sx q[3];
rz(2.485399) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.31768826) q[0];
sx q[0];
rz(-2.3264139) q[0];
sx q[0];
rz(1.748388) q[0];
rz(-2.4488917) q[1];
sx q[1];
rz(-1.0857948) q[1];
sx q[1];
rz(1.0903953) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4596953) q[0];
sx q[0];
rz(-2.5715264) q[0];
sx q[0];
rz(1.1718114) q[0];
x q[1];
rz(2.6165025) q[2];
sx q[2];
rz(-1.0193362) q[2];
sx q[2];
rz(-0.22071018) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-2.322256) q[1];
sx q[1];
rz(-2.6988479) q[1];
sx q[1];
rz(1.2135452) q[1];
rz(-pi) q[2];
rz(-0.15655915) q[3];
sx q[3];
rz(-1.0362175) q[3];
sx q[3];
rz(1.8434356) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(2.580487) q[2];
sx q[2];
rz(-1.1539536) q[2];
sx q[2];
rz(0.53606501) q[2];
rz(-2.8481893) q[3];
sx q[3];
rz(-1.2111726) q[3];
sx q[3];
rz(-0.48040473) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2368161) q[0];
sx q[0];
rz(-0.95229709) q[0];
sx q[0];
rz(-0.099040898) q[0];
rz(1.0320484) q[1];
sx q[1];
rz(-0.85056225) q[1];
sx q[1];
rz(1.7031857) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8598762) q[0];
sx q[0];
rz(-1.3163573) q[0];
sx q[0];
rz(-0.86812015) q[0];
rz(-pi) q[1];
rz(-1.4969669) q[2];
sx q[2];
rz(-1.2880688) q[2];
sx q[2];
rz(-1.7746995) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-0.26353961) q[1];
sx q[1];
rz(-0.78070736) q[1];
sx q[1];
rz(-1.6040463) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.766633) q[3];
sx q[3];
rz(-2.8228033) q[3];
sx q[3];
rz(2.6998883) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-1.9859887) q[2];
sx q[2];
rz(-0.5039379) q[2];
sx q[2];
rz(-0.8209374) q[2];
rz(1.692903) q[3];
sx q[3];
rz(-2.4235348) q[3];
sx q[3];
rz(-1.4887571) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.89707017) q[0];
sx q[0];
rz(-2.2269766) q[0];
sx q[0];
rz(2.6389417) q[0];
rz(-1.2333168) q[1];
sx q[1];
rz(-0.93524593) q[1];
sx q[1];
rz(0.9800235) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.21979688) q[0];
sx q[0];
rz(-1.3338519) q[0];
sx q[0];
rz(1.3150929) q[0];
rz(-0.29666846) q[2];
sx q[2];
rz(-0.7985332) q[2];
sx q[2];
rz(0.980033) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(2.7777694) q[1];
sx q[1];
rz(-2.2968493) q[1];
sx q[1];
rz(-1.2593624) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.3408889) q[3];
sx q[3];
rz(-2.3977444) q[3];
sx q[3];
rz(-0.4901674) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(0.22214733) q[2];
sx q[2];
rz(-1.8541226) q[2];
sx q[2];
rz(-1.7670828) q[2];
rz(-1.3162656) q[3];
sx q[3];
rz(-2.2856183) q[3];
sx q[3];
rz(-0.98178664) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.114349) q[0];
sx q[0];
rz(-2.7482432) q[0];
sx q[0];
rz(-2.223176) q[0];
rz(2.9367327) q[1];
sx q[1];
rz(-1.058895) q[1];
sx q[1];
rz(-1.6580261) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.61648864) q[0];
sx q[0];
rz(-1.8601928) q[0];
sx q[0];
rz(0.15286907) q[0];
x q[1];
rz(2.3307033) q[2];
sx q[2];
rz(-1.6605018) q[2];
sx q[2];
rz(0.16505884) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-3.0682448) q[1];
sx q[1];
rz(-0.7052583) q[1];
sx q[1];
rz(1.7068935) q[1];
rz(-0.056413944) q[3];
sx q[3];
rz(-0.71995633) q[3];
sx q[3];
rz(-0.82701339) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(2.1398937) q[2];
sx q[2];
rz(-2.6504982) q[2];
sx q[2];
rz(-1.8074544) q[2];
rz(0.88207465) q[3];
sx q[3];
rz(-0.69552723) q[3];
sx q[3];
rz(-1.2923366) q[3];
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
x q[0];
rz(pi/2) q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0598711) q[0];
sx q[0];
rz(-2.7099755) q[0];
sx q[0];
rz(3.1234142) q[0];
rz(-0.40953088) q[1];
sx q[1];
rz(-1.4298226) q[1];
sx q[1];
rz(3.0407564) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4977328) q[0];
sx q[0];
rz(-2.5063305) q[0];
sx q[0];
rz(-0.058202581) q[0];
rz(-pi) q[1];
rz(0.42068215) q[2];
sx q[2];
rz(-1.5245588) q[2];
sx q[2];
rz(1.3476404) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.455115) q[1];
sx q[1];
rz(-1.2077189) q[1];
sx q[1];
rz(-2.8836714) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.12091095) q[3];
sx q[3];
rz(-1.8163346) q[3];
sx q[3];
rz(-2.3335378) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-0.33710256) q[2];
sx q[2];
rz(-3.0323995) q[2];
sx q[2];
rz(-1.8936554) q[2];
rz(-1.2248056) q[3];
sx q[3];
rz(-0.8453415) q[3];
sx q[3];
rz(3.0758408) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9058022) q[0];
sx q[0];
rz(-1.6658655) q[0];
sx q[0];
rz(0.32050785) q[0];
rz(2.7505752) q[1];
sx q[1];
rz(-2.0000439) q[1];
sx q[1];
rz(-1.5919707) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.81330106) q[0];
sx q[0];
rz(-1.3718954) q[0];
sx q[0];
rz(1.3764253) q[0];
x q[1];
rz(-2.8587927) q[2];
sx q[2];
rz(-0.53164266) q[2];
sx q[2];
rz(-0.34257364) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(2.446881) q[1];
sx q[1];
rz(-1.2031735) q[1];
sx q[1];
rz(0.093631677) q[1];
rz(-0.016716047) q[3];
sx q[3];
rz(-2.7191945) q[3];
sx q[3];
rz(1.0836386) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(1.090342) q[2];
sx q[2];
rz(-1.0498468) q[2];
sx q[2];
rz(2.7412097) q[2];
rz(-2.9054437) q[3];
sx q[3];
rz(-0.103424) q[3];
sx q[3];
rz(2.1307438) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.878933) q[0];
sx q[0];
rz(-2.6317609) q[0];
sx q[0];
rz(1.2608933) q[0];
rz(-0.48514584) q[1];
sx q[1];
rz(-0.79882516) q[1];
sx q[1];
rz(-0.47732236) q[1];
rz(-1.7791228) q[2];
sx q[2];
rz(-1.4478417) q[2];
sx q[2];
rz(1.7461591) q[2];
rz(-1.7988206) q[3];
sx q[3];
rz(-1.6327471) q[3];
sx q[3];
rz(-3.0625797) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
