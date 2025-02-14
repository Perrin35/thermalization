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
rz(-1.7464632) q[0];
sx q[0];
rz(-0.75463086) q[0];
sx q[0];
rz(-2.057743) q[0];
rz(2.2639182) q[1];
sx q[1];
rz(4.35507) q[1];
sx q[1];
rz(9.9014643) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3221368) q[0];
sx q[0];
rz(-1.1753723) q[0];
sx q[0];
rz(-0.28015341) q[0];
rz(-pi) q[1];
rz(2.0332912) q[2];
sx q[2];
rz(-1.1574189) q[2];
sx q[2];
rz(0.55962901) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-1.4922766) q[1];
sx q[1];
rz(-1.5879803) q[1];
sx q[1];
rz(1.8284998) q[1];
rz(-pi) q[2];
rz(-1.3506557) q[3];
sx q[3];
rz(-1.1486037) q[3];
sx q[3];
rz(3.0289502) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-3.1276663) q[2];
sx q[2];
rz(-1.4907336) q[2];
sx q[2];
rz(0.15065436) q[2];
rz(-1.6023747) q[3];
sx q[3];
rz(-2.7645002) q[3];
sx q[3];
rz(-0.25130513) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(pi/2) q[3];
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
rz(-0.056331228) q[0];
sx q[0];
rz(-1.769861) q[0];
sx q[0];
rz(2.7714609) q[0];
rz(2.6689957) q[1];
sx q[1];
rz(-0.31817803) q[1];
sx q[1];
rz(0.95742375) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3963378) q[0];
sx q[0];
rz(-0.50512132) q[0];
sx q[0];
rz(1.4046679) q[0];
x q[1];
rz(-1.1195807) q[2];
sx q[2];
rz(-1.7282082) q[2];
sx q[2];
rz(1.6119286) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-0.70499252) q[1];
sx q[1];
rz(-2.3967522) q[1];
sx q[1];
rz(0.1257433) q[1];
x q[2];
rz(-1.9022835) q[3];
sx q[3];
rz(-0.9210862) q[3];
sx q[3];
rz(0.099516221) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-0.0051673278) q[2];
sx q[2];
rz(-2.397126) q[2];
sx q[2];
rz(1.0235419) q[2];
rz(-1.3905904) q[3];
sx q[3];
rz(-1.3955045) q[3];
sx q[3];
rz(1.8243779) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
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
rz(0.35109529) q[0];
sx q[0];
rz(-2.1320765) q[0];
sx q[0];
rz(0.56280953) q[0];
rz(-1.6273181) q[1];
sx q[1];
rz(-1.4311675) q[1];
sx q[1];
rz(-0.079158457) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1437627) q[0];
sx q[0];
rz(-2.061686) q[0];
sx q[0];
rz(1.2116371) q[0];
rz(-pi) q[1];
rz(1.6116033) q[2];
sx q[2];
rz(-1.0804324) q[2];
sx q[2];
rz(0.23863579) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(2.1677163) q[1];
sx q[1];
rz(-1.9922246) q[1];
sx q[1];
rz(-2.3371479) q[1];
rz(-pi) q[2];
x q[2];
rz(2.2658562) q[3];
sx q[3];
rz(-0.99627021) q[3];
sx q[3];
rz(0.72573001) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(1.9074273) q[2];
sx q[2];
rz(-1.5831466) q[2];
sx q[2];
rz(-1.3321715) q[2];
rz(2.7857156) q[3];
sx q[3];
rz(-1.9477113) q[3];
sx q[3];
rz(-2.1607384) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6818162) q[0];
sx q[0];
rz(-1.1695319) q[0];
sx q[0];
rz(-2.0303149) q[0];
rz(-2.2214644) q[1];
sx q[1];
rz(-1.5524813) q[1];
sx q[1];
rz(0.2535325) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.64719114) q[0];
sx q[0];
rz(-1.8888942) q[0];
sx q[0];
rz(2.6107016) q[0];
rz(-pi) q[1];
rz(-2.9811834) q[2];
sx q[2];
rz(-1.4575851) q[2];
sx q[2];
rz(-2.8594174) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-0.82099709) q[1];
sx q[1];
rz(-1.7554469) q[1];
sx q[1];
rz(0.76652938) q[1];
rz(-2.8315968) q[3];
sx q[3];
rz(-1.0862203) q[3];
sx q[3];
rz(-2.692401) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(0.8186875) q[2];
sx q[2];
rz(-0.56537586) q[2];
sx q[2];
rz(1.4595855) q[2];
rz(-2.5679576) q[3];
sx q[3];
rz(-1.2021659) q[3];
sx q[3];
rz(-3.0972163) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6673073) q[0];
sx q[0];
rz(-2.0720338) q[0];
sx q[0];
rz(1.3473508) q[0];
rz(-0.59066331) q[1];
sx q[1];
rz(-1.9213516) q[1];
sx q[1];
rz(1.4264872) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.41353696) q[0];
sx q[0];
rz(-1.1794834) q[0];
sx q[0];
rz(1.9584459) q[0];
rz(-pi) q[1];
rz(-0.64138647) q[2];
sx q[2];
rz(-1.8872627) q[2];
sx q[2];
rz(-3.0015025) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(3.0795143) q[1];
sx q[1];
rz(-0.2582363) q[1];
sx q[1];
rz(-2.9596427) q[1];
rz(-pi) q[2];
x q[2];
rz(1.5065881) q[3];
sx q[3];
rz(-1.9252464) q[3];
sx q[3];
rz(-1.3469157) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(2.9023989) q[2];
sx q[2];
rz(-2.0812483) q[2];
sx q[2];
rz(1.9690465) q[2];
rz(-0.086056195) q[3];
sx q[3];
rz(-0.9525758) q[3];
sx q[3];
rz(-0.15611592) q[3];
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
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0758783) q[0];
sx q[0];
rz(-1.5944163) q[0];
sx q[0];
rz(-0.86268798) q[0];
rz(0.24738303) q[1];
sx q[1];
rz(-1.8377973) q[1];
sx q[1];
rz(0.47201306) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2675572) q[0];
sx q[0];
rz(-2.3568601) q[0];
sx q[0];
rz(-2.6536056) q[0];
x q[1];
rz(-1.3026779) q[2];
sx q[2];
rz(-2.6202218) q[2];
sx q[2];
rz(0.15586317) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.0942877) q[1];
sx q[1];
rz(-1.3940812) q[1];
sx q[1];
rz(-0.30728886) q[1];
x q[2];
rz(1.8627164) q[3];
sx q[3];
rz(-2.4987767) q[3];
sx q[3];
rz(1.1228648) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(0.90999675) q[2];
sx q[2];
rz(-2.1746077) q[2];
sx q[2];
rz(1.7426573) q[2];
rz(1.4553962) q[3];
sx q[3];
rz(-2.3536436) q[3];
sx q[3];
rz(-0.56070352) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6550734) q[0];
sx q[0];
rz(-1.1460679) q[0];
sx q[0];
rz(2.6229677) q[0];
rz(-0.71867603) q[1];
sx q[1];
rz(-0.51281723) q[1];
sx q[1];
rz(1.3291043) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.14940093) q[0];
sx q[0];
rz(-1.2103545) q[0];
sx q[0];
rz(-2.0130231) q[0];
rz(-pi) q[1];
rz(-1.7460004) q[2];
sx q[2];
rz(-0.28131286) q[2];
sx q[2];
rz(-3.0228721) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(0.79256637) q[1];
sx q[1];
rz(-1.6036803) q[1];
sx q[1];
rz(-2.3057641) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.5409327) q[3];
sx q[3];
rz(-2.5549508) q[3];
sx q[3];
rz(-1.3636953) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.79080498) q[2];
sx q[2];
rz(-0.69956508) q[2];
sx q[2];
rz(0.038185509) q[2];
rz(-2.3112678) q[3];
sx q[3];
rz(-1.7540951) q[3];
sx q[3];
rz(-3.0939046) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4627948) q[0];
sx q[0];
rz(-2.9252453) q[0];
sx q[0];
rz(1.3702673) q[0];
rz(-0.38617745) q[1];
sx q[1];
rz(-1.736707) q[1];
sx q[1];
rz(-2.4716299) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7012335) q[0];
sx q[0];
rz(-1.0770849) q[0];
sx q[0];
rz(2.6026653) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.4248895) q[2];
sx q[2];
rz(-1.6316902) q[2];
sx q[2];
rz(-0.74541192) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(1.252287) q[1];
sx q[1];
rz(-2.5355967) q[1];
sx q[1];
rz(1.3713981) q[1];
rz(-pi) q[2];
x q[2];
rz(3.1042323) q[3];
sx q[3];
rz(-1.3589763) q[3];
sx q[3];
rz(-2.7461014) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(2.2204444) q[2];
sx q[2];
rz(-0.28382742) q[2];
sx q[2];
rz(-1.05668) q[2];
rz(-1.4165261) q[3];
sx q[3];
rz(-1.2456015) q[3];
sx q[3];
rz(-1.2904803) q[3];
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
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0678299) q[0];
sx q[0];
rz(-2.1011598) q[0];
sx q[0];
rz(2.1347866) q[0];
rz(-2.6489068) q[1];
sx q[1];
rz(-1.410781) q[1];
sx q[1];
rz(2.3115092) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0706129) q[0];
sx q[0];
rz(-1.604053) q[0];
sx q[0];
rz(-0.279056) q[0];
rz(0.30980457) q[2];
sx q[2];
rz(-2.5199157) q[2];
sx q[2];
rz(-0.72337389) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(0.54129564) q[1];
sx q[1];
rz(-1.7233371) q[1];
sx q[1];
rz(-0.18220724) q[1];
x q[2];
rz(1.5437627) q[3];
sx q[3];
rz(-0.44315954) q[3];
sx q[3];
rz(-1.5704607) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(0.81862187) q[2];
sx q[2];
rz(-1.3016204) q[2];
sx q[2];
rz(2.0446365) q[2];
rz(-1.8638301) q[3];
sx q[3];
rz(-0.68456972) q[3];
sx q[3];
rz(-2.4975615) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.85035664) q[0];
sx q[0];
rz(-2.0391897) q[0];
sx q[0];
rz(-0.37503234) q[0];
rz(0.11861435) q[1];
sx q[1];
rz(-1.7053441) q[1];
sx q[1];
rz(2.2172701) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3548022) q[0];
sx q[0];
rz(-1.2806007) q[0];
sx q[0];
rz(-0.59450392) q[0];
rz(2.5293328) q[2];
sx q[2];
rz(-2.3697402) q[2];
sx q[2];
rz(1.5764232) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(1.3872747) q[1];
sx q[1];
rz(-1.5994834) q[1];
sx q[1];
rz(-0.36092511) q[1];
rz(-pi) q[2];
rz(-0.64226867) q[3];
sx q[3];
rz(-0.54107252) q[3];
sx q[3];
rz(-2.9716345) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-2.2529926) q[2];
sx q[2];
rz(-2.1414976) q[2];
sx q[2];
rz(0.89317733) q[2];
rz(2.0231953) q[3];
sx q[3];
rz(-1.018254) q[3];
sx q[3];
rz(0.65677381) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7669582) q[0];
sx q[0];
rz(-1.081291) q[0];
sx q[0];
rz(-1.3670856) q[0];
rz(0.18192667) q[1];
sx q[1];
rz(-0.55768273) q[1];
sx q[1];
rz(2.2348977) q[1];
rz(1.8858269) q[2];
sx q[2];
rz(-1.4834822) q[2];
sx q[2];
rz(-0.33374141) q[2];
rz(-0.51707324) q[3];
sx q[3];
rz(-1.2096249) q[3];
sx q[3];
rz(0.67740868) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
