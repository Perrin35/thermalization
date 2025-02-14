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
rz(1.9407152) q[0];
sx q[0];
rz(-2.241029) q[0];
sx q[0];
rz(0.23298921) q[0];
rz(1.6147344) q[1];
sx q[1];
rz(-1.072071) q[1];
sx q[1];
rz(-1.1195247) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6699864) q[0];
sx q[0];
rz(-1.3363839) q[0];
sx q[0];
rz(0.013201518) q[0];
x q[1];
rz(-0.018466516) q[2];
sx q[2];
rz(-0.56150836) q[2];
sx q[2];
rz(0.36040053) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-0.2762766) q[1];
sx q[1];
rz(-2.7807693) q[1];
sx q[1];
rz(1.3555384) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.60910881) q[3];
sx q[3];
rz(-1.2909596) q[3];
sx q[3];
rz(2.5697054) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(2.4423674) q[2];
sx q[2];
rz(-2.0825601) q[2];
sx q[2];
rz(-3.0287058) q[2];
rz(-0.26116192) q[3];
sx q[3];
rz(-1.7966725) q[3];
sx q[3];
rz(1.2507218) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.79171044) q[0];
sx q[0];
rz(-0.078502027) q[0];
sx q[0];
rz(0.054340266) q[0];
rz(0.21866523) q[1];
sx q[1];
rz(-1.4726787) q[1];
sx q[1];
rz(-2.7770619) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2549627) q[0];
sx q[0];
rz(-0.97575649) q[0];
sx q[0];
rz(0.9592077) q[0];
rz(-pi) q[1];
x q[1];
rz(2.0979375) q[2];
sx q[2];
rz(-2.4424565) q[2];
sx q[2];
rz(1.3329891) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-1.5082701) q[1];
sx q[1];
rz(-2.5145217) q[1];
sx q[1];
rz(3.0116664) q[1];
x q[2];
rz(-2.0516615) q[3];
sx q[3];
rz(-0.37853795) q[3];
sx q[3];
rz(2.7906281) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-2.9008122) q[2];
sx q[2];
rz(-1.6996926) q[2];
sx q[2];
rz(1.1085294) q[2];
rz(2.945914) q[3];
sx q[3];
rz(-1.207573) q[3];
sx q[3];
rz(1.4669363) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.3942669) q[0];
sx q[0];
rz(-2.1013923) q[0];
sx q[0];
rz(1.0362097) q[0];
rz(-2.5569083) q[1];
sx q[1];
rz(-1.6744637) q[1];
sx q[1];
rz(-0.55589693) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.79093874) q[0];
sx q[0];
rz(-2.2430111) q[0];
sx q[0];
rz(2.1124798) q[0];
x q[1];
rz(1.118752) q[2];
sx q[2];
rz(-3.0556745) q[2];
sx q[2];
rz(-0.8762067) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.30196721) q[1];
sx q[1];
rz(-2.4263067) q[1];
sx q[1];
rz(-1.3705181) q[1];
rz(1.317765) q[3];
sx q[3];
rz(-1.1339427) q[3];
sx q[3];
rz(-1.4405314) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(1.3802152) q[2];
sx q[2];
rz(-0.20647241) q[2];
sx q[2];
rz(-1.192344) q[2];
rz(-0.10382593) q[3];
sx q[3];
rz(-1.1622279) q[3];
sx q[3];
rz(-1.1473568) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.1059859) q[0];
sx q[0];
rz(-0.62653956) q[0];
sx q[0];
rz(-1.71738) q[0];
rz(0.72319889) q[1];
sx q[1];
rz(-2.9056748) q[1];
sx q[1];
rz(0.22044388) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0766692) q[0];
sx q[0];
rz(-0.60416302) q[0];
sx q[0];
rz(-0.88978617) q[0];
rz(2.1313558) q[2];
sx q[2];
rz(-1.3187485) q[2];
sx q[2];
rz(0.28365669) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-0.38762906) q[1];
sx q[1];
rz(-2.7668608) q[1];
sx q[1];
rz(1.8829569) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.8414874) q[3];
sx q[3];
rz(-0.53181767) q[3];
sx q[3];
rz(1.8790661) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(0.38349884) q[2];
sx q[2];
rz(-1.3145964) q[2];
sx q[2];
rz(-2.626075) q[2];
rz(3.0564485) q[3];
sx q[3];
rz(-2.4862423) q[3];
sx q[3];
rz(0.83318025) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.67007095) q[0];
sx q[0];
rz(-2.9586198) q[0];
sx q[0];
rz(-0.66437379) q[0];
rz(2.6145256) q[1];
sx q[1];
rz(-2.0030237) q[1];
sx q[1];
rz(-1.1767496) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7648286) q[0];
sx q[0];
rz(-1.485678) q[0];
sx q[0];
rz(-0.028193817) q[0];
rz(0.88567733) q[2];
sx q[2];
rz(-1.9373477) q[2];
sx q[2];
rz(-1.4057856) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(2.5866111) q[1];
sx q[1];
rz(-2.9025295) q[1];
sx q[1];
rz(-3.0596517) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.9620911) q[3];
sx q[3];
rz(-2.13509) q[3];
sx q[3];
rz(1.5998942) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(1.1145757) q[2];
sx q[2];
rz(-0.21275529) q[2];
sx q[2];
rz(2.2034755) q[2];
rz(1.045687) q[3];
sx q[3];
rz(-0.18200471) q[3];
sx q[3];
rz(-0.81030455) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi/2) q[3];
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
rz(-0.806458) q[0];
sx q[0];
rz(-1.9430176) q[0];
sx q[0];
rz(-1.8916116) q[0];
rz(-0.20420034) q[1];
sx q[1];
rz(-0.67664346) q[1];
sx q[1];
rz(2.2883889) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.51914224) q[0];
sx q[0];
rz(-1.2850437) q[0];
sx q[0];
rz(0.35730548) q[0];
rz(-pi) q[1];
rz(2.3647652) q[2];
sx q[2];
rz(-2.3558801) q[2];
sx q[2];
rz(-1.7091027) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-2.0029162) q[1];
sx q[1];
rz(-2.2670548) q[1];
sx q[1];
rz(0.55714861) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.12746396) q[3];
sx q[3];
rz(-1.2317622) q[3];
sx q[3];
rz(1.6515428) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.044067232) q[2];
sx q[2];
rz(-1.7996457) q[2];
sx q[2];
rz(-1.1996783) q[2];
rz(2.1357644) q[3];
sx q[3];
rz(-0.89322105) q[3];
sx q[3];
rz(-0.83542663) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.47578874) q[0];
sx q[0];
rz(-2.8732193) q[0];
sx q[0];
rz(-0.98261181) q[0];
rz(1.8334552) q[1];
sx q[1];
rz(-1.2160701) q[1];
sx q[1];
rz(1.4391724) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.84837259) q[0];
sx q[0];
rz(-1.2407899) q[0];
sx q[0];
rz(2.7820935) q[0];
rz(-pi) q[1];
x q[1];
rz(1.0487532) q[2];
sx q[2];
rz(-1.5035275) q[2];
sx q[2];
rz(-1.0053588) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(0.01165302) q[1];
sx q[1];
rz(-2.2352152) q[1];
sx q[1];
rz(2.7205977) q[1];
rz(-0.47726008) q[3];
sx q[3];
rz(-1.4443732) q[3];
sx q[3];
rz(1.3341158) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(2.9787489) q[2];
sx q[2];
rz(-2.2480201) q[2];
sx q[2];
rz(2.5992375) q[2];
rz(-0.98313037) q[3];
sx q[3];
rz(-1.977908) q[3];
sx q[3];
rz(0.94016176) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5777609) q[0];
sx q[0];
rz(-1.4262119) q[0];
sx q[0];
rz(2.3610709) q[0];
rz(-1.966194) q[1];
sx q[1];
rz(-0.54882097) q[1];
sx q[1];
rz(-1.0343879) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7808409) q[0];
sx q[0];
rz(-2.7324317) q[0];
sx q[0];
rz(-0.43171127) q[0];
rz(2.0883191) q[2];
sx q[2];
rz(-0.1642326) q[2];
sx q[2];
rz(0.17271481) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.24496291) q[1];
sx q[1];
rz(-1.1707998) q[1];
sx q[1];
rz(-1.6952312) q[1];
rz(2.5938321) q[3];
sx q[3];
rz(-1.3410853) q[3];
sx q[3];
rz(-1.821777) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-0.15022755) q[2];
sx q[2];
rz(-1.7346202) q[2];
sx q[2];
rz(0.11037174) q[2];
rz(-1.8433833) q[3];
sx q[3];
rz(-1.7896264) q[3];
sx q[3];
rz(2.8492294) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.89123911) q[0];
sx q[0];
rz(-1.013101) q[0];
sx q[0];
rz(-1.891834) q[0];
rz(3.0505782) q[1];
sx q[1];
rz(-0.72151557) q[1];
sx q[1];
rz(0.7116085) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2219951) q[0];
sx q[0];
rz(-1.5775497) q[0];
sx q[0];
rz(-2.3632977) q[0];
x q[1];
rz(0.83699147) q[2];
sx q[2];
rz(-0.92304936) q[2];
sx q[2];
rz(-0.68488065) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(1.8597128) q[1];
sx q[1];
rz(-1.2285474) q[1];
sx q[1];
rz(2.6564391) q[1];
x q[2];
rz(-1.2290088) q[3];
sx q[3];
rz(-0.84484378) q[3];
sx q[3];
rz(-1.0085229) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-2.8083501) q[2];
sx q[2];
rz(-1.4834206) q[2];
sx q[2];
rz(1.9688155) q[2];
rz(1.6138389) q[3];
sx q[3];
rz(-2.0634191) q[3];
sx q[3];
rz(0.095002256) q[3];
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
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.67977366) q[0];
sx q[0];
rz(-0.12400308) q[0];
sx q[0];
rz(-1.5661731) q[0];
rz(2.7836986) q[1];
sx q[1];
rz(-1.0708258) q[1];
sx q[1];
rz(-2.2391052) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1514777) q[0];
sx q[0];
rz(-0.95989908) q[0];
sx q[0];
rz(0.31661242) q[0];
rz(-pi) q[1];
x q[1];
rz(1.5487368) q[2];
sx q[2];
rz(-2.214693) q[2];
sx q[2];
rz(-0.48286706) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-1.5229228) q[1];
sx q[1];
rz(-0.82397193) q[1];
sx q[1];
rz(2.0342779) q[1];
x q[2];
rz(-1.7125413) q[3];
sx q[3];
rz(-0.55484164) q[3];
sx q[3];
rz(-2.6920163) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-2.0365399) q[2];
sx q[2];
rz(-0.5883216) q[2];
sx q[2];
rz(-1.6416637) q[2];
rz(-2.6203652) q[3];
sx q[3];
rz(-0.14385496) q[3];
sx q[3];
rz(-2.1874793) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.70984107) q[0];
sx q[0];
rz(-2.3713645) q[0];
sx q[0];
rz(1.4159528) q[0];
rz(-3.1406828) q[1];
sx q[1];
rz(-1.4716499) q[1];
sx q[1];
rz(1.6843527) q[1];
rz(1.7224689) q[2];
sx q[2];
rz(-0.57874741) q[2];
sx q[2];
rz(-0.49606965) q[2];
rz(2.4215563) q[3];
sx q[3];
rz(-2.2514718) q[3];
sx q[3];
rz(2.9801647) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
