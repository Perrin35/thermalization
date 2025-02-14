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
rz(-0.44218818) q[0];
sx q[0];
rz(-1.1531885) q[0];
sx q[0];
rz(3.0129504) q[0];
rz(-4.6075912) q[1];
sx q[1];
rz(1.0990376) q[1];
sx q[1];
rz(8.3429835) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2483467) q[0];
sx q[0];
rz(-1.5061646) q[0];
sx q[0];
rz(-1.6016225) q[0];
x q[1];
rz(-1.9456909) q[2];
sx q[2];
rz(-1.8842976) q[2];
sx q[2];
rz(1.9913265) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-1.9777367) q[1];
sx q[1];
rz(-2.9675936) q[1];
sx q[1];
rz(-1.15088) q[1];
x q[2];
rz(-3.0469037) q[3];
sx q[3];
rz(-1.296954) q[3];
sx q[3];
rz(1.3665592) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-3.0470524) q[2];
sx q[2];
rz(-0.69539842) q[2];
sx q[2];
rz(-0.73235861) q[2];
rz(0.098946027) q[3];
sx q[3];
rz(-2.1061335) q[3];
sx q[3];
rz(-1.8100479) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.089379646) q[0];
sx q[0];
rz(-2.3591924) q[0];
sx q[0];
rz(0.041444929) q[0];
rz(-2.4502358) q[1];
sx q[1];
rz(-1.3203011) q[1];
sx q[1];
rz(2.2533805) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.81538686) q[0];
sx q[0];
rz(-1.0727709) q[0];
sx q[0];
rz(2.2603358) q[0];
rz(-pi) q[1];
rz(0.58063796) q[2];
sx q[2];
rz(-2.3742445) q[2];
sx q[2];
rz(1.5754171) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.2728558) q[1];
sx q[1];
rz(-1.5275729) q[1];
sx q[1];
rz(-2.8985913) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.51155297) q[3];
sx q[3];
rz(-1.8163876) q[3];
sx q[3];
rz(-3.0893774) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-1.9899675) q[2];
sx q[2];
rz(-1.7513559) q[2];
sx q[2];
rz(-1.7612696) q[2];
rz(-0.049792854) q[3];
sx q[3];
rz(-2.4536665) q[3];
sx q[3];
rz(-0.5602347) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3675073) q[0];
sx q[0];
rz(-2.9893576) q[0];
sx q[0];
rz(1.2138858) q[0];
rz(-1.8269352) q[1];
sx q[1];
rz(-2.1036802) q[1];
sx q[1];
rz(-1.5710057) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9938068) q[0];
sx q[0];
rz(-1.8166313) q[0];
sx q[0];
rz(-2.7509806) q[0];
rz(-pi) q[1];
x q[1];
rz(0.41061398) q[2];
sx q[2];
rz(-1.7143692) q[2];
sx q[2];
rz(-0.5856572) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(1.5017058) q[1];
sx q[1];
rz(-2.9738131) q[1];
sx q[1];
rz(1.3659992) q[1];
rz(-pi) q[2];
rz(2.5475471) q[3];
sx q[3];
rz(-1.0540773) q[3];
sx q[3];
rz(-1.3498211) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-2.7870002) q[2];
sx q[2];
rz(-2.3921693) q[2];
sx q[2];
rz(-2.836239) q[2];
rz(-2.2954156) q[3];
sx q[3];
rz(-1.2140707) q[3];
sx q[3];
rz(-2.2727216) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2884035) q[0];
sx q[0];
rz(-2.1222293) q[0];
sx q[0];
rz(0.97002059) q[0];
rz(-2.8963529) q[1];
sx q[1];
rz(-1.0673362) q[1];
sx q[1];
rz(1.6397569) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.95432276) q[0];
sx q[0];
rz(-1.6980972) q[0];
sx q[0];
rz(2.9856647) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.5108893) q[2];
sx q[2];
rz(-0.28894385) q[2];
sx q[2];
rz(-2.9693791) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(1.3093331) q[1];
sx q[1];
rz(-2.8974779) q[1];
sx q[1];
rz(-2.4629618) q[1];
x q[2];
rz(1.8957696) q[3];
sx q[3];
rz(-0.87034908) q[3];
sx q[3];
rz(-2.9440299) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(2.5331427) q[2];
sx q[2];
rz(-1.5773982) q[2];
sx q[2];
rz(1.8428141) q[2];
rz(-0.53457824) q[3];
sx q[3];
rz(-0.60408533) q[3];
sx q[3];
rz(0.66876137) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.088148549) q[0];
sx q[0];
rz(-1.7867333) q[0];
sx q[0];
rz(-2.163929) q[0];
rz(-2.4560302) q[1];
sx q[1];
rz(-2.3473163) q[1];
sx q[1];
rz(0.96644863) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9622877) q[0];
sx q[0];
rz(-2.0317269) q[0];
sx q[0];
rz(-0.75051744) q[0];
x q[1];
rz(-0.62018779) q[2];
sx q[2];
rz(-0.66588565) q[2];
sx q[2];
rz(-2.263139) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(0.19627608) q[1];
sx q[1];
rz(-1.4320201) q[1];
sx q[1];
rz(-0.86244418) q[1];
rz(-pi) q[2];
rz(2.590476) q[3];
sx q[3];
rz(-1.2663906) q[3];
sx q[3];
rz(-1.9929043) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.55296772) q[2];
sx q[2];
rz(-1.3543465) q[2];
sx q[2];
rz(0.28263131) q[2];
rz(2.1770554) q[3];
sx q[3];
rz(-1.5805565) q[3];
sx q[3];
rz(2.7717822) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
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
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.10837567) q[0];
sx q[0];
rz(-0.41154698) q[0];
sx q[0];
rz(-2.0630398) q[0];
rz(-2.318553) q[1];
sx q[1];
rz(-1.1230725) q[1];
sx q[1];
rz(-0.80026921) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5259994) q[0];
sx q[0];
rz(-1.5999927) q[0];
sx q[0];
rz(0.008511624) q[0];
rz(-pi) q[1];
rz(-2.8326464) q[2];
sx q[2];
rz(-1.2703203) q[2];
sx q[2];
rz(0.38602513) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-2.2507402) q[1];
sx q[1];
rz(-1.0167334) q[1];
sx q[1];
rz(-2.3046369) q[1];
rz(-pi) q[2];
x q[2];
rz(2.5979302) q[3];
sx q[3];
rz(-2.398536) q[3];
sx q[3];
rz(2.0280968) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.1396973) q[2];
sx q[2];
rz(-1.2418094) q[2];
sx q[2];
rz(-0.099420698) q[2];
rz(0.56441489) q[3];
sx q[3];
rz(-1.5921009) q[3];
sx q[3];
rz(-2.2216589) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.99892202) q[0];
sx q[0];
rz(-0.26695928) q[0];
sx q[0];
rz(-0.024918407) q[0];
rz(-2.0492367) q[1];
sx q[1];
rz(-1.9905636) q[1];
sx q[1];
rz(0.54164642) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3643862) q[0];
sx q[0];
rz(-1.9310675) q[0];
sx q[0];
rz(2.2132753) q[0];
rz(-1.0850026) q[2];
sx q[2];
rz(-0.62720229) q[2];
sx q[2];
rz(-1.021334) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(0.4723963) q[1];
sx q[1];
rz(-1.7150567) q[1];
sx q[1];
rz(-1.902182) q[1];
rz(-pi) q[2];
x q[2];
rz(1.8239519) q[3];
sx q[3];
rz(-2.5196155) q[3];
sx q[3];
rz(3.1282116) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-0.8267374) q[2];
sx q[2];
rz(-1.9172226) q[2];
sx q[2];
rz(1.1401736) q[2];
rz(-1.7013811) q[3];
sx q[3];
rz(-0.39339104) q[3];
sx q[3];
rz(-2.9668729) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2575271) q[0];
sx q[0];
rz(-1.9767569) q[0];
sx q[0];
rz(2.2648532) q[0];
rz(-2.9675617) q[1];
sx q[1];
rz(-2.048025) q[1];
sx q[1];
rz(1.3678975) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.68326658) q[0];
sx q[0];
rz(-1.84671) q[0];
sx q[0];
rz(1.0066731) q[0];
rz(-pi) q[1];
rz(2.9139745) q[2];
sx q[2];
rz(-1.2151067) q[2];
sx q[2];
rz(1.5700036) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(1.9570043) q[1];
sx q[1];
rz(-2.4589029) q[1];
sx q[1];
rz(0.95889389) q[1];
rz(-pi) q[2];
x q[2];
rz(0.075966751) q[3];
sx q[3];
rz(-1.6923762) q[3];
sx q[3];
rz(1.0978804) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-1.6690663) q[2];
sx q[2];
rz(-1.3269227) q[2];
sx q[2];
rz(0.18829045) q[2];
rz(-1.3663728) q[3];
sx q[3];
rz(-1.3683189) q[3];
sx q[3];
rz(-1.9395456) q[3];
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
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.61757225) q[0];
sx q[0];
rz(-0.13372788) q[0];
sx q[0];
rz(-2.4753841) q[0];
rz(1.1353525) q[1];
sx q[1];
rz(-0.98061457) q[1];
sx q[1];
rz(1.1599783) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.85765094) q[0];
sx q[0];
rz(-0.91383024) q[0];
sx q[0];
rz(-0.82842555) q[0];
rz(-pi) q[1];
rz(0.23187821) q[2];
sx q[2];
rz(-2.3740951) q[2];
sx q[2];
rz(-0.19814834) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-3.0409909) q[1];
sx q[1];
rz(-1.5631277) q[1];
sx q[1];
rz(2.1867382) q[1];
rz(1.0053535) q[3];
sx q[3];
rz(-1.3015038) q[3];
sx q[3];
rz(-2.5776742) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.74405115) q[2];
sx q[2];
rz(-0.56637374) q[2];
sx q[2];
rz(-2.9446824) q[2];
rz(-2.2714553) q[3];
sx q[3];
rz(-2.0700442) q[3];
sx q[3];
rz(-1.8740694) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5043735) q[0];
sx q[0];
rz(-2.0688031) q[0];
sx q[0];
rz(-3.1043501) q[0];
rz(1.0875018) q[1];
sx q[1];
rz(-2.0934413) q[1];
sx q[1];
rz(-0.16669272) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8213538) q[0];
sx q[0];
rz(-1.3036089) q[0];
sx q[0];
rz(0.70810476) q[0];
x q[1];
rz(0.48515452) q[2];
sx q[2];
rz(-1.8388766) q[2];
sx q[2];
rz(-1.9659276) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.31754759) q[1];
sx q[1];
rz(-0.30774857) q[1];
sx q[1];
rz(-2.8466562) q[1];
rz(0.44346614) q[3];
sx q[3];
rz(-2.9812668) q[3];
sx q[3];
rz(3.0976618) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.4166261) q[2];
sx q[2];
rz(-0.61369696) q[2];
sx q[2];
rz(-1.7566768) q[2];
rz(-2.7409605) q[3];
sx q[3];
rz(-1.8568361) q[3];
sx q[3];
rz(1.0111151) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8818125) q[0];
sx q[0];
rz(-2.0406944) q[0];
sx q[0];
rz(2.5115321) q[0];
rz(0.13314816) q[1];
sx q[1];
rz(-1.9825736) q[1];
sx q[1];
rz(2.0864743) q[1];
rz(1.7693949) q[2];
sx q[2];
rz(-0.19459859) q[2];
sx q[2];
rz(2.9698402) q[2];
rz(0.12689982) q[3];
sx q[3];
rz(-0.77346934) q[3];
sx q[3];
rz(1.0897549) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
