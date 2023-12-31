OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-1.3553319) q[0];
sx q[0];
rz(0.064602764) q[0];
sx q[0];
rz(6.3048007) q[0];
rz(2.1463483) q[1];
sx q[1];
rz(-1.8145476) q[1];
sx q[1];
rz(-1.8099161) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.74237139) q[0];
sx q[0];
rz(-2.135015) q[0];
sx q[0];
rz(-2.7447634) q[0];
rz(-pi) q[1];
rz(-0.10279074) q[2];
sx q[2];
rz(-1.79709) q[2];
sx q[2];
rz(3.0068827) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-1.761258) q[1];
sx q[1];
rz(-0.77612703) q[1];
sx q[1];
rz(2.9556729) q[1];
x q[2];
rz(-2.8033923) q[3];
sx q[3];
rz(-0.97465289) q[3];
sx q[3];
rz(2.3232834) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-2.8109479) q[2];
sx q[2];
rz(-1.4928879) q[2];
sx q[2];
rz(2.5374106) q[2];
rz(1.0243246) q[3];
sx q[3];
rz(-1.9842792) q[3];
sx q[3];
rz(-2.0143051) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
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
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8169096) q[0];
sx q[0];
rz(-0.01318251) q[0];
sx q[0];
rz(-1.0634134) q[0];
rz(-2.2564607) q[1];
sx q[1];
rz(-1.5567895) q[1];
sx q[1];
rz(0.0016454776) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.20365276) q[0];
sx q[0];
rz(-1.5814039) q[0];
sx q[0];
rz(1.6158576) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.33866377) q[2];
sx q[2];
rz(-0.67064697) q[2];
sx q[2];
rz(1.2884969) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(2.8811223) q[1];
sx q[1];
rz(-0.80658856) q[1];
sx q[1];
rz(-1.9116198) q[1];
x q[2];
rz(0.9886338) q[3];
sx q[3];
rz(-2.8391264) q[3];
sx q[3];
rz(-0.31103381) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.985618) q[2];
sx q[2];
rz(-1.6321471) q[2];
sx q[2];
rz(-2.4334811) q[2];
rz(-2.0478785) q[3];
sx q[3];
rz(-1.8023068) q[3];
sx q[3];
rz(-0.99350199) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.920632) q[0];
sx q[0];
rz(-1.7376124) q[0];
sx q[0];
rz(2.3143342) q[0];
rz(-3.1365085) q[1];
sx q[1];
rz(-1.2132443) q[1];
sx q[1];
rz(2.0522096) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.1339061) q[0];
sx q[0];
rz(-1.0431164) q[0];
sx q[0];
rz(-0.49467996) q[0];
x q[1];
rz(1.1459648) q[2];
sx q[2];
rz(-1.0152738) q[2];
sx q[2];
rz(-0.59031634) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(0.032420302) q[1];
sx q[1];
rz(-1.65616) q[1];
sx q[1];
rz(1.1847772) q[1];
rz(-pi) q[2];
x q[2];
rz(2.5855519) q[3];
sx q[3];
rz(-2.2696113) q[3];
sx q[3];
rz(1.3211105) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-2.0777145) q[2];
sx q[2];
rz(-2.2641247) q[2];
sx q[2];
rz(-2.2568978) q[2];
rz(1.1832773) q[3];
sx q[3];
rz(-1.3180472) q[3];
sx q[3];
rz(-1.4900835) q[3];
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
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5692212) q[0];
sx q[0];
rz(-0.5643934) q[0];
sx q[0];
rz(-0.72682056) q[0];
rz(-0.80365333) q[1];
sx q[1];
rz(-1.0875965) q[1];
sx q[1];
rz(2.7817536) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7774178) q[0];
sx q[0];
rz(-2.4628277) q[0];
sx q[0];
rz(-1.2749519) q[0];
rz(1.6845409) q[2];
sx q[2];
rz(-0.92278381) q[2];
sx q[2];
rz(-1.5930454) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-2.847536) q[1];
sx q[1];
rz(-1.9901853) q[1];
sx q[1];
rz(-1.9860752) q[1];
rz(-pi) q[2];
x q[2];
rz(0.8538586) q[3];
sx q[3];
rz(-2.6918908) q[3];
sx q[3];
rz(-2.2886697) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(0.56746733) q[2];
sx q[2];
rz(-1.987477) q[2];
sx q[2];
rz(0.7652258) q[2];
rz(-2.3848173) q[3];
sx q[3];
rz(-0.59195834) q[3];
sx q[3];
rz(-0.24100196) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9969479) q[0];
sx q[0];
rz(-0.47676555) q[0];
sx q[0];
rz(-1.7255406) q[0];
rz(0.36711806) q[1];
sx q[1];
rz(-1.7558302) q[1];
sx q[1];
rz(-1.0353154) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2238335) q[0];
sx q[0];
rz(-0.98063722) q[0];
sx q[0];
rz(-2.607164) q[0];
rz(-pi) q[1];
rz(0.073280235) q[2];
sx q[2];
rz(-0.49256941) q[2];
sx q[2];
rz(-0.57893334) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-0.1086515) q[1];
sx q[1];
rz(-1.6254394) q[1];
sx q[1];
rz(-3.0461237) q[1];
rz(-pi) q[2];
rz(-1.7807998) q[3];
sx q[3];
rz(-1.8404507) q[3];
sx q[3];
rz(0.3013914) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(1.3570024) q[2];
sx q[2];
rz(-1.4214397) q[2];
sx q[2];
rz(-0.76888293) q[2];
rz(2.8055577) q[3];
sx q[3];
rz(-0.78032812) q[3];
sx q[3];
rz(-1.4679573) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1795905) q[0];
sx q[0];
rz(-0.63350326) q[0];
sx q[0];
rz(1.1451716) q[0];
rz(1.1046474) q[1];
sx q[1];
rz(-1.9112174) q[1];
sx q[1];
rz(2.9343658) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4743054) q[0];
sx q[0];
rz(-1.7100088) q[0];
sx q[0];
rz(2.5179203) q[0];
rz(-pi) q[1];
rz(1.1057165) q[2];
sx q[2];
rz(-1.0208703) q[2];
sx q[2];
rz(-1.20649) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(2.0761557) q[1];
sx q[1];
rz(-1.3216615) q[1];
sx q[1];
rz(0.25030987) q[1];
rz(1.3979785) q[3];
sx q[3];
rz(-0.76068766) q[3];
sx q[3];
rz(-2.2067604) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.8032288) q[2];
sx q[2];
rz(-0.88023606) q[2];
sx q[2];
rz(-0.72171372) q[2];
rz(-1.2747814) q[3];
sx q[3];
rz(-1.5694247) q[3];
sx q[3];
rz(2.8815564) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8969144) q[0];
sx q[0];
rz(-1.5667863) q[0];
sx q[0];
rz(-2.4095643) q[0];
rz(0.0094982068) q[1];
sx q[1];
rz(-2.5932725) q[1];
sx q[1];
rz(2.8498555) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0193298) q[0];
sx q[0];
rz(-2.7195192) q[0];
sx q[0];
rz(-0.12515573) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.2675915) q[2];
sx q[2];
rz(-1.469194) q[2];
sx q[2];
rz(-2.6333957) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(1.8818672) q[1];
sx q[1];
rz(-1.7181267) q[1];
sx q[1];
rz(-0.30585652) q[1];
rz(1.8802059) q[3];
sx q[3];
rz(-2.072812) q[3];
sx q[3];
rz(1.3537784) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(2.8774524) q[2];
sx q[2];
rz(-1.143127) q[2];
sx q[2];
rz(-0.83731246) q[2];
rz(1.9705747) q[3];
sx q[3];
rz(-1.5191017) q[3];
sx q[3];
rz(3.0795735) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
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
rz(-1.0432805) q[0];
sx q[0];
rz(-0.68798143) q[0];
sx q[0];
rz(-2.4777381) q[0];
rz(-3.035416) q[1];
sx q[1];
rz(-0.60634923) q[1];
sx q[1];
rz(0.95867872) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7784568) q[0];
sx q[0];
rz(-2.083337) q[0];
sx q[0];
rz(-2.1666359) q[0];
rz(-pi) q[1];
rz(1.1218698) q[2];
sx q[2];
rz(-0.58510963) q[2];
sx q[2];
rz(-1.4657071) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-2.9740323) q[1];
sx q[1];
rz(-1.4640199) q[1];
sx q[1];
rz(2.0978931) q[1];
rz(-pi) q[2];
rz(-2.9577191) q[3];
sx q[3];
rz(-0.5920147) q[3];
sx q[3];
rz(1.5541935) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(2.133193) q[2];
sx q[2];
rz(-1.2352751) q[2];
sx q[2];
rz(-2.2224902) q[2];
rz(1.7637926) q[3];
sx q[3];
rz(-1.2104687) q[3];
sx q[3];
rz(1.1834043) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8063426) q[0];
sx q[0];
rz(-0.85334539) q[0];
sx q[0];
rz(2.7046955) q[0];
rz(-2.4412952) q[1];
sx q[1];
rz(-1.4236139) q[1];
sx q[1];
rz(2.6760496) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.31637329) q[0];
sx q[0];
rz(-1.6557949) q[0];
sx q[0];
rz(-2.4337016) q[0];
rz(-pi) q[1];
rz(-1.1586458) q[2];
sx q[2];
rz(-1.4530384) q[2];
sx q[2];
rz(-2.462537) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(2.000396) q[1];
sx q[1];
rz(-1.2277514) q[1];
sx q[1];
rz(-2.6617906) q[1];
rz(1.9302619) q[3];
sx q[3];
rz(-1.8065479) q[3];
sx q[3];
rz(2.204493) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.7312701) q[2];
sx q[2];
rz(-2.6022544) q[2];
sx q[2];
rz(-1.643606) q[2];
rz(0.20478976) q[3];
sx q[3];
rz(-1.0054761) q[3];
sx q[3];
rz(-0.27206102) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4964504) q[0];
sx q[0];
rz(-2.5800939) q[0];
sx q[0];
rz(1.1219332) q[0];
rz(-2.3902068) q[1];
sx q[1];
rz(-1.6151927) q[1];
sx q[1];
rz(0.5823935) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.84401417) q[0];
sx q[0];
rz(-1.3676924) q[0];
sx q[0];
rz(-0.88589478) q[0];
rz(-0.15162823) q[2];
sx q[2];
rz(-1.0851589) q[2];
sx q[2];
rz(-3.1079353) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.460234) q[1];
sx q[1];
rz(-1.6885307) q[1];
sx q[1];
rz(1.387499) q[1];
rz(3.0155229) q[3];
sx q[3];
rz(-2.9226916) q[3];
sx q[3];
rz(-0.88760469) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.0104684) q[2];
sx q[2];
rz(-0.98103395) q[2];
sx q[2];
rz(-0.76254145) q[2];
rz(-1.4108346) q[3];
sx q[3];
rz(-0.91791955) q[3];
sx q[3];
rz(-1.5677174) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0762155) q[0];
sx q[0];
rz(-2.1552754) q[0];
sx q[0];
rz(-1.7393204) q[0];
rz(-1.3394042) q[1];
sx q[1];
rz(-1.6911472) q[1];
sx q[1];
rz(-1.4858248) q[1];
rz(0.81417685) q[2];
sx q[2];
rz(-1.3635175) q[2];
sx q[2];
rz(-1.2843532) q[2];
rz(0.81621697) q[3];
sx q[3];
rz(-0.64707884) q[3];
sx q[3];
rz(-0.94221471) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
