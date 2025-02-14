OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-2.7944827) q[0];
sx q[0];
rz(-1.4248983) q[0];
sx q[0];
rz(0.27118924) q[0];
rz(-2.7746692) q[1];
sx q[1];
rz(-0.75777268) q[1];
sx q[1];
rz(-1.4581207) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9957033) q[0];
sx q[0];
rz(-1.5738537) q[0];
sx q[0];
rz(1.5637335) q[0];
rz(1.0881937) q[2];
sx q[2];
rz(-1.6763684) q[2];
sx q[2];
rz(0.93868773) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-0.62195146) q[1];
sx q[1];
rz(-0.84712815) q[1];
sx q[1];
rz(-2.3569941) q[1];
rz(-pi) q[2];
x q[2];
rz(-3.0786282) q[3];
sx q[3];
rz(-0.28423542) q[3];
sx q[3];
rz(0.96719826) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(0.014341982) q[2];
sx q[2];
rz(-0.6659826) q[2];
sx q[2];
rz(-1.2464397) q[2];
rz(-1.0154826) q[3];
sx q[3];
rz(-0.4963488) q[3];
sx q[3];
rz(2.523876) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.604973) q[0];
sx q[0];
rz(-0.07285694) q[0];
sx q[0];
rz(2.4541722) q[0];
rz(0.61126002) q[1];
sx q[1];
rz(-1.6300853) q[1];
sx q[1];
rz(-0.90868178) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.1686104) q[0];
sx q[0];
rz(-1.6618722) q[0];
sx q[0];
rz(-0.9206207) q[0];
rz(0.15678997) q[2];
sx q[2];
rz(-2.1794381) q[2];
sx q[2];
rz(-1.551924) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-0.37088739) q[1];
sx q[1];
rz(-1.0418278) q[1];
sx q[1];
rz(-2.4338399) q[1];
x q[2];
rz(0.55697316) q[3];
sx q[3];
rz(-1.4629629) q[3];
sx q[3];
rz(1.7086017) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.0574135) q[2];
sx q[2];
rz(-0.97687352) q[2];
sx q[2];
rz(1.9834391) q[2];
rz(0.17364764) q[3];
sx q[3];
rz(-1.4519139) q[3];
sx q[3];
rz(0.68918532) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3678906) q[0];
sx q[0];
rz(-2.9359718) q[0];
sx q[0];
rz(-0.50262991) q[0];
rz(1.8058808) q[1];
sx q[1];
rz(-3.0323961) q[1];
sx q[1];
rz(1.4044382) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3356467) q[0];
sx q[0];
rz(-2.5281161) q[0];
sx q[0];
rz(-2.4015266) q[0];
rz(-pi) q[1];
rz(2.7680306) q[2];
sx q[2];
rz(-1.3907916) q[2];
sx q[2];
rz(1.9587868) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(2.6413414) q[1];
sx q[1];
rz(-1.6585357) q[1];
sx q[1];
rz(0.32890202) q[1];
rz(-2.5504774) q[3];
sx q[3];
rz(-1.1151259) q[3];
sx q[3];
rz(-2.299813) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-1.3001083) q[2];
sx q[2];
rz(-0.69977641) q[2];
sx q[2];
rz(2.4669199) q[2];
rz(-2.789433) q[3];
sx q[3];
rz(-1.9904741) q[3];
sx q[3];
rz(-1.6262416) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.801055) q[0];
sx q[0];
rz(-0.041943701) q[0];
sx q[0];
rz(1.1658143) q[0];
rz(-0.55533987) q[1];
sx q[1];
rz(-0.97709877) q[1];
sx q[1];
rz(-1.4110483) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4395907) q[0];
sx q[0];
rz(-1.5541663) q[0];
sx q[0];
rz(2.8921333) q[0];
rz(-pi) q[1];
x q[1];
rz(2.0981815) q[2];
sx q[2];
rz(-1.4031271) q[2];
sx q[2];
rz(2.4737918) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-0.41563161) q[1];
sx q[1];
rz(-1.2228194) q[1];
sx q[1];
rz(-1.0370952) q[1];
rz(-pi) q[2];
rz(-0.23616236) q[3];
sx q[3];
rz(-2.3104226) q[3];
sx q[3];
rz(2.1787852) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(2.1022776) q[2];
sx q[2];
rz(-0.5216051) q[2];
sx q[2];
rz(-2.6585141) q[2];
rz(-2.3704902) q[3];
sx q[3];
rz(-1.5604138) q[3];
sx q[3];
rz(1.34351) q[3];
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
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1191206) q[0];
sx q[0];
rz(-2.6799057) q[0];
sx q[0];
rz(-0.22892496) q[0];
rz(-1.6772894) q[1];
sx q[1];
rz(-2.5720282) q[1];
sx q[1];
rz(-0.48008188) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.10413607) q[0];
sx q[0];
rz(-2.7711282) q[0];
sx q[0];
rz(-2.3568826) q[0];
x q[1];
rz(-1.3111021) q[2];
sx q[2];
rz(-1.9995481) q[2];
sx q[2];
rz(1.4285806) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-0.11949524) q[1];
sx q[1];
rz(-1.198631) q[1];
sx q[1];
rz(1.5161901) q[1];
rz(-2.9579406) q[3];
sx q[3];
rz(-0.95377659) q[3];
sx q[3];
rz(2.4434095) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(3.0698801) q[2];
sx q[2];
rz(-2.0166848) q[2];
sx q[2];
rz(1.2360392) q[2];
rz(-1.49336) q[3];
sx q[3];
rz(-0.83087102) q[3];
sx q[3];
rz(1.4504455) q[3];
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
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.013783197) q[0];
sx q[0];
rz(-2.4007128) q[0];
sx q[0];
rz(0.86454779) q[0];
rz(-2.3132482) q[1];
sx q[1];
rz(-1.8048077) q[1];
sx q[1];
rz(0.72992647) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2695539) q[0];
sx q[0];
rz(-2.3667496) q[0];
sx q[0];
rz(-0.86685769) q[0];
rz(-pi) q[1];
rz(0.66940825) q[2];
sx q[2];
rz(-1.727549) q[2];
sx q[2];
rz(2.0956958) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(1.1600274) q[1];
sx q[1];
rz(-1.3277819) q[1];
sx q[1];
rz(1.0286938) q[1];
rz(-pi) q[2];
x q[2];
rz(-3.0804068) q[3];
sx q[3];
rz(-2.4172999) q[3];
sx q[3];
rz(-0.35288399) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(0.56968969) q[2];
sx q[2];
rz(-0.14293417) q[2];
sx q[2];
rz(-1.8611056) q[2];
rz(1.5754383) q[3];
sx q[3];
rz(-2.1011293) q[3];
sx q[3];
rz(-0.86281002) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.60523072) q[0];
sx q[0];
rz(-1.0030168) q[0];
sx q[0];
rz(-0.59148106) q[0];
rz(-2.483967) q[1];
sx q[1];
rz(-1.769519) q[1];
sx q[1];
rz(-3.0234911) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.5786835) q[0];
sx q[0];
rz(-1.9171035) q[0];
sx q[0];
rz(-1.3265394) q[0];
rz(-0.12848358) q[2];
sx q[2];
rz(-2.621935) q[2];
sx q[2];
rz(-2.3574587) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.98886314) q[1];
sx q[1];
rz(-2.0901457) q[1];
sx q[1];
rz(1.6822271) q[1];
rz(-0.16987993) q[3];
sx q[3];
rz(-1.3515444) q[3];
sx q[3];
rz(-0.77283054) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-2.9401231) q[2];
sx q[2];
rz(-2.757299) q[2];
sx q[2];
rz(1.8817792) q[2];
rz(-1.2093557) q[3];
sx q[3];
rz(-1.7989379) q[3];
sx q[3];
rz(-0.847009) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.92846337) q[0];
sx q[0];
rz(-2.6452112) q[0];
sx q[0];
rz(-1.553836) q[0];
rz(1.3262879) q[1];
sx q[1];
rz(-0.98618788) q[1];
sx q[1];
rz(2.3415668) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.532142) q[0];
sx q[0];
rz(-2.6649619) q[0];
sx q[0];
rz(0.6269639) q[0];
rz(-pi) q[1];
x q[1];
rz(0.50234143) q[2];
sx q[2];
rz(-2.3141344) q[2];
sx q[2];
rz(2.7934144) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-3.0985089) q[1];
sx q[1];
rz(-1.0810163) q[1];
sx q[1];
rz(-2.5950762) q[1];
rz(-pi) q[2];
rz(-2.7553237) q[3];
sx q[3];
rz(-1.0479386) q[3];
sx q[3];
rz(2.8012365) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(2.942261) q[2];
sx q[2];
rz(-2.3776725) q[2];
sx q[2];
rz(-1.2660816) q[2];
rz(2.074504) q[3];
sx q[3];
rz(-1.4184003) q[3];
sx q[3];
rz(-0.057723109) q[3];
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
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.98457321) q[0];
sx q[0];
rz(-0.83681256) q[0];
sx q[0];
rz(2.9465604) q[0];
rz(0.86206478) q[1];
sx q[1];
rz(-1.74086) q[1];
sx q[1];
rz(-0.21924266) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.040960366) q[0];
sx q[0];
rz(-2.1233243) q[0];
sx q[0];
rz(-0.29968963) q[0];
rz(-1.1226038) q[2];
sx q[2];
rz(-1.7328129) q[2];
sx q[2];
rz(-2.3336099) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(0.91867204) q[1];
sx q[1];
rz(-1.5377475) q[1];
sx q[1];
rz(2.6085684) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.3211649) q[3];
sx q[3];
rz(-1.9747989) q[3];
sx q[3];
rz(-1.8694316) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(1.4054883) q[2];
sx q[2];
rz(-0.89724237) q[2];
sx q[2];
rz(1.7264977) q[2];
rz(-1.803558) q[3];
sx q[3];
rz(-1.3058563) q[3];
sx q[3];
rz(0.068917902) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
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
rz(-0.59704429) q[0];
sx q[0];
rz(-2.2291849) q[0];
sx q[0];
rz(1.9507116) q[0];
rz(1.4471588) q[1];
sx q[1];
rz(-1.8444599) q[1];
sx q[1];
rz(1.7097293) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.53729355) q[0];
sx q[0];
rz(-1.7794637) q[0];
sx q[0];
rz(0.66402004) q[0];
x q[1];
rz(-2.5578021) q[2];
sx q[2];
rz(-2.4707831) q[2];
sx q[2];
rz(-1.5354615) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.398795) q[1];
sx q[1];
rz(-0.94091641) q[1];
sx q[1];
rz(2.6349318) q[1];
rz(-pi) q[2];
x q[2];
rz(1.6064243) q[3];
sx q[3];
rz(-2.4089353) q[3];
sx q[3];
rz(-3.0095524) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-0.65149629) q[2];
sx q[2];
rz(-2.8054674) q[2];
sx q[2];
rz(2.2641342) q[2];
rz(-1.6089926) q[3];
sx q[3];
rz(-1.420615) q[3];
sx q[3];
rz(2.1144313) q[3];
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
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0090016) q[0];
sx q[0];
rz(-1.4132211) q[0];
sx q[0];
rz(2.6268517) q[0];
rz(-2.4280836) q[1];
sx q[1];
rz(-2.3206354) q[1];
sx q[1];
rz(-1.6101507) q[1];
rz(-0.64856957) q[2];
sx q[2];
rz(-2.6852457) q[2];
sx q[2];
rz(-1.1732875) q[2];
rz(0.71047495) q[3];
sx q[3];
rz(-2.8024891) q[3];
sx q[3];
rz(2.9935318) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
