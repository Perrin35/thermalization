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
rz(-3.0769899) q[0];
sx q[0];
rz(3.1199772) q[0];
rz(-0.99524438) q[1];
sx q[1];
rz(-1.3270451) q[1];
sx q[1];
rz(1.8099161) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0639122) q[0];
sx q[0];
rz(-2.4644116) q[0];
sx q[0];
rz(-2.1190686) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.7982499) q[2];
sx q[2];
rz(-1.4706352) q[2];
sx q[2];
rz(-1.7286466) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.3803346) q[1];
sx q[1];
rz(-2.3654656) q[1];
sx q[1];
rz(0.18591979) q[1];
rz(0.94727912) q[3];
sx q[3];
rz(-1.2926658) q[3];
sx q[3];
rz(0.55752414) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(2.8109479) q[2];
sx q[2];
rz(-1.6487048) q[2];
sx q[2];
rz(2.5374106) q[2];
rz(-2.1172681) q[3];
sx q[3];
rz(-1.1573135) q[3];
sx q[3];
rz(-1.1272875) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.3246831) q[0];
sx q[0];
rz(-3.1284101) q[0];
sx q[0];
rz(1.0634134) q[0];
rz(-2.2564607) q[1];
sx q[1];
rz(-1.5848031) q[1];
sx q[1];
rz(-0.0016454776) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.20365276) q[0];
sx q[0];
rz(-1.5814039) q[0];
sx q[0];
rz(1.6158576) q[0];
rz(-pi) q[1];
rz(2.4992141) q[2];
sx q[2];
rz(-1.3628236) q[2];
sx q[2];
rz(3.1285398) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(1.5509847) q[1];
sx q[1];
rz(-1.814517) q[1];
sx q[1];
rz(2.3477712) q[1];
rz(-pi) q[2];
x q[2];
rz(1.8257636) q[3];
sx q[3];
rz(-1.4062738) q[3];
sx q[3];
rz(0.69873519) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(1.985618) q[2];
sx q[2];
rz(-1.5094455) q[2];
sx q[2];
rz(-0.70811159) q[2];
rz(-2.0478785) q[3];
sx q[3];
rz(-1.3392859) q[3];
sx q[3];
rz(-2.1480907) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.920632) q[0];
sx q[0];
rz(-1.4039803) q[0];
sx q[0];
rz(2.3143342) q[0];
rz(3.1365085) q[1];
sx q[1];
rz(-1.9283483) q[1];
sx q[1];
rz(-1.089383) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.0076865772) q[0];
sx q[0];
rz(-1.0431164) q[0];
sx q[0];
rz(2.6469127) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.59805869) q[2];
sx q[2];
rz(-1.9285678) q[2];
sx q[2];
rz(1.2146815) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-1.3315862) q[1];
sx q[1];
rz(-0.3948822) q[1];
sx q[1];
rz(-1.3473131) q[1];
x q[2];
rz(0.79077625) q[3];
sx q[3];
rz(-1.1547935) q[3];
sx q[3];
rz(2.5115867) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-2.0777145) q[2];
sx q[2];
rz(-2.2641247) q[2];
sx q[2];
rz(0.88469488) q[2];
rz(-1.9583154) q[3];
sx q[3];
rz(-1.8235455) q[3];
sx q[3];
rz(1.4900835) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5723715) q[0];
sx q[0];
rz(-2.5771993) q[0];
sx q[0];
rz(-0.72682056) q[0];
rz(0.80365333) q[1];
sx q[1];
rz(-2.0539961) q[1];
sx q[1];
rz(-0.35983905) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.737405) q[0];
sx q[0];
rz(-0.92659896) q[0];
sx q[0];
rz(-2.9106211) q[0];
rz(-pi) q[1];
rz(1.6845409) q[2];
sx q[2];
rz(-0.92278381) q[2];
sx q[2];
rz(-1.5930454) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.6872014) q[1];
sx q[1];
rz(-1.9481716) q[1];
sx q[1];
rz(0.45339938) q[1];
rz(0.8538586) q[3];
sx q[3];
rz(-0.44970185) q[3];
sx q[3];
rz(2.2886697) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-2.5741253) q[2];
sx q[2];
rz(-1.987477) q[2];
sx q[2];
rz(0.7652258) q[2];
rz(-0.75677538) q[3];
sx q[3];
rz(-2.5496343) q[3];
sx q[3];
rz(-0.24100196) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1446447) q[0];
sx q[0];
rz(-0.47676555) q[0];
sx q[0];
rz(1.416052) q[0];
rz(-0.36711806) q[1];
sx q[1];
rz(-1.3857625) q[1];
sx q[1];
rz(-1.0353154) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.1128164) q[0];
sx q[0];
rz(-2.007764) q[0];
sx q[0];
rz(-0.90941888) q[0];
x q[1];
rz(1.6100699) q[2];
sx q[2];
rz(-1.0796667) q[2];
sx q[2];
rz(0.66205762) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.1086515) q[1];
sx q[1];
rz(-1.5161533) q[1];
sx q[1];
rz(3.0461237) q[1];
x q[2];
rz(1.7807998) q[3];
sx q[3];
rz(-1.3011419) q[3];
sx q[3];
rz(0.3013914) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(1.3570024) q[2];
sx q[2];
rz(-1.4214397) q[2];
sx q[2];
rz(-0.76888293) q[2];
rz(-0.33603493) q[3];
sx q[3];
rz(-0.78032812) q[3];
sx q[3];
rz(1.6736354) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.96200213) q[0];
sx q[0];
rz(-0.63350326) q[0];
sx q[0];
rz(1.996421) q[0];
rz(-1.1046474) q[1];
sx q[1];
rz(-1.9112174) q[1];
sx q[1];
rz(0.2072269) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4743054) q[0];
sx q[0];
rz(-1.7100088) q[0];
sx q[0];
rz(2.5179203) q[0];
rz(-pi) q[1];
x q[1];
rz(0.60116641) q[2];
sx q[2];
rz(-1.9631557) q[2];
sx q[2];
rz(2.5208134) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(1.065437) q[1];
sx q[1];
rz(-1.3216615) q[1];
sx q[1];
rz(-0.25030987) q[1];
rz(1.3979785) q[3];
sx q[3];
rz(-0.76068766) q[3];
sx q[3];
rz(0.93483227) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.3383639) q[2];
sx q[2];
rz(-2.2613566) q[2];
sx q[2];
rz(2.4198789) q[2];
rz(1.8668113) q[3];
sx q[3];
rz(-1.5721679) q[3];
sx q[3];
rz(0.26003626) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8969144) q[0];
sx q[0];
rz(-1.5748064) q[0];
sx q[0];
rz(0.73202837) q[0];
rz(-3.1320944) q[1];
sx q[1];
rz(-0.54832012) q[1];
sx q[1];
rz(0.2917372) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7043982) q[0];
sx q[0];
rz(-1.5196374) q[0];
sx q[0];
rz(-0.41914661) q[0];
x q[1];
rz(-1.2675915) q[2];
sx q[2];
rz(-1.6723987) q[2];
sx q[2];
rz(-0.50819699) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(1.8818672) q[1];
sx q[1];
rz(-1.7181267) q[1];
sx q[1];
rz(-0.30585652) q[1];
rz(-pi) q[2];
rz(-1.8802059) q[3];
sx q[3];
rz(-2.072812) q[3];
sx q[3];
rz(1.7878143) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.26414028) q[2];
sx q[2];
rz(-1.9984657) q[2];
sx q[2];
rz(-2.3042802) q[2];
rz(-1.9705747) q[3];
sx q[3];
rz(-1.6224909) q[3];
sx q[3];
rz(3.0795735) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
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
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0432805) q[0];
sx q[0];
rz(-2.4536112) q[0];
sx q[0];
rz(-0.6638546) q[0];
rz(-3.035416) q[1];
sx q[1];
rz(-2.5352434) q[1];
sx q[1];
rz(2.1829139) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7230969) q[0];
sx q[0];
rz(-2.3765058) q[0];
sx q[0];
rz(0.78406783) q[0];
rz(-pi) q[1];
rz(1.1218698) q[2];
sx q[2];
rz(-0.58510963) q[2];
sx q[2];
rz(1.6758855) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.9740323) q[1];
sx q[1];
rz(-1.6775727) q[1];
sx q[1];
rz(2.0978931) q[1];
rz(2.9577191) q[3];
sx q[3];
rz(-2.549578) q[3];
sx q[3];
rz(-1.5873991) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-1.0083996) q[2];
sx q[2];
rz(-1.9063176) q[2];
sx q[2];
rz(2.2224902) q[2];
rz(1.7637926) q[3];
sx q[3];
rz(-1.2104687) q[3];
sx q[3];
rz(1.1834043) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8063426) q[0];
sx q[0];
rz(-0.85334539) q[0];
sx q[0];
rz(-0.43689716) q[0];
rz(0.70029744) q[1];
sx q[1];
rz(-1.4236139) q[1];
sx q[1];
rz(-0.46554309) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3532928) q[0];
sx q[0];
rz(-0.71209891) q[0];
sx q[0];
rz(-0.13029356) q[0];
x q[1];
rz(1.9829468) q[2];
sx q[2];
rz(-1.4530384) q[2];
sx q[2];
rz(0.6790557) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-2.000396) q[1];
sx q[1];
rz(-1.2277514) q[1];
sx q[1];
rz(2.6617906) q[1];
rz(-1.9302619) q[3];
sx q[3];
rz(-1.3350447) q[3];
sx q[3];
rz(2.204493) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.41032252) q[2];
sx q[2];
rz(-0.53933829) q[2];
sx q[2];
rz(1.643606) q[2];
rz(0.20478976) q[3];
sx q[3];
rz(-2.1361165) q[3];
sx q[3];
rz(-2.8695316) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4964504) q[0];
sx q[0];
rz(-0.56149879) q[0];
sx q[0];
rz(-1.1219332) q[0];
rz(-2.3902068) q[1];
sx q[1];
rz(-1.6151927) q[1];
sx q[1];
rz(0.5823935) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.56349194) q[0];
sx q[0];
rz(-0.90257114) q[0];
sx q[0];
rz(-2.8816954) q[0];
x q[1];
rz(-2.0612129) q[2];
sx q[2];
rz(-1.4368125) q[2];
sx q[2];
rz(1.4659363) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(3.0092587) q[1];
sx q[1];
rz(-1.3887822) q[1];
sx q[1];
rz(-0.11972129) q[1];
rz(-pi) q[2];
x q[2];
rz(2.9243745) q[3];
sx q[3];
rz(-1.5434885) q[3];
sx q[3];
rz(2.581493) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-2.1311243) q[2];
sx q[2];
rz(-0.98103395) q[2];
sx q[2];
rz(0.76254145) q[2];
rz(-1.7307581) q[3];
sx q[3];
rz(-0.91791955) q[3];
sx q[3];
rz(1.5677174) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0653771) q[0];
sx q[0];
rz(-2.1552754) q[0];
sx q[0];
rz(-1.7393204) q[0];
rz(-1.8021884) q[1];
sx q[1];
rz(-1.4504455) q[1];
sx q[1];
rz(1.6557678) q[1];
rz(1.273524) q[2];
sx q[2];
rz(-0.77902972) q[2];
sx q[2];
rz(-3.0697889) q[2];
rz(-1.0675666) q[3];
sx q[3];
rz(-1.9964841) q[3];
sx q[3];
rz(-1.869429) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
