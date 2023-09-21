OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-2.6887309) q[0];
sx q[0];
rz(-2.9714669) q[0];
sx q[0];
rz(-2.3556019) q[0];
rz(0.6056447) q[1];
sx q[1];
rz(3.7925386) q[1];
sx q[1];
rz(8.7906919) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2274322) q[0];
sx q[0];
rz(-2.3254546) q[0];
sx q[0];
rz(0.039365191) q[0];
x q[1];
rz(1.6526821) q[2];
sx q[2];
rz(-2.4856644) q[2];
sx q[2];
rz(-1.3002849) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(0.69971985) q[1];
sx q[1];
rz(-2.3708214) q[1];
sx q[1];
rz(-2.2295879) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.4001873) q[3];
sx q[3];
rz(-1.549198) q[3];
sx q[3];
rz(-0.17822972) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-1.9156076) q[2];
sx q[2];
rz(-0.51708022) q[2];
sx q[2];
rz(-1.263164) q[2];
rz(1.2849215) q[3];
sx q[3];
rz(-1.4611171) q[3];
sx q[3];
rz(-2.8485956) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
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
rz(0.15853515) q[0];
sx q[0];
rz(-1.2936658) q[0];
sx q[0];
rz(-2.7040226) q[0];
rz(0.63105398) q[1];
sx q[1];
rz(-0.32866207) q[1];
sx q[1];
rz(0.22110573) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3723345) q[0];
sx q[0];
rz(-0.49494574) q[0];
sx q[0];
rz(-1.9702205) q[0];
rz(-pi) q[1];
rz(1.2220076) q[2];
sx q[2];
rz(-1.4663327) q[2];
sx q[2];
rz(-1.2806569) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-1.4202538) q[1];
sx q[1];
rz(-1.4495279) q[1];
sx q[1];
rz(-0.46055693) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.50237327) q[3];
sx q[3];
rz(-1.4851735) q[3];
sx q[3];
rz(1.3115713) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-0.21800403) q[2];
sx q[2];
rz(-1.4322586) q[2];
sx q[2];
rz(0.58829266) q[2];
rz(-2.6925987) q[3];
sx q[3];
rz(-2.7136927) q[3];
sx q[3];
rz(2.0935521) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.56851971) q[0];
sx q[0];
rz(-0.097671106) q[0];
sx q[0];
rz(1.0004689) q[0];
rz(0.72552848) q[1];
sx q[1];
rz(-1.0745445) q[1];
sx q[1];
rz(-0.75769889) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.66540668) q[0];
sx q[0];
rz(-1.2520257) q[0];
sx q[0];
rz(-0.55253367) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.4865815) q[2];
sx q[2];
rz(-0.25946028) q[2];
sx q[2];
rz(-0.26817817) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(2.0211027) q[1];
sx q[1];
rz(-2.0334525) q[1];
sx q[1];
rz(0.91336577) q[1];
rz(1.4593616) q[3];
sx q[3];
rz(-2.206344) q[3];
sx q[3];
rz(0.089126822) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.9329325) q[2];
sx q[2];
rz(-2.9186086) q[2];
sx q[2];
rz(-2.3051252) q[2];
rz(1.6992735) q[3];
sx q[3];
rz(-1.6675555) q[3];
sx q[3];
rz(2.9377655) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[3];
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
rz(0.02012415) q[0];
sx q[0];
rz(-3.0614873) q[0];
sx q[0];
rz(-2.6304723) q[0];
rz(0.077443667) q[1];
sx q[1];
rz(-0.42235342) q[1];
sx q[1];
rz(-1.8130594) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7610361) q[0];
sx q[0];
rz(-1.3576926) q[0];
sx q[0];
rz(-0.15844945) q[0];
rz(-pi) q[1];
x q[1];
rz(0.51283299) q[2];
sx q[2];
rz(-2.3094258) q[2];
sx q[2];
rz(1.2733449) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.37175825) q[1];
sx q[1];
rz(-1.8375988) q[1];
sx q[1];
rz(-3.1410602) q[1];
x q[2];
rz(1.7926746) q[3];
sx q[3];
rz(-1.4639877) q[3];
sx q[3];
rz(-2.2282003) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(0.093734309) q[2];
sx q[2];
rz(-1.9231984) q[2];
sx q[2];
rz(1.1070586) q[2];
rz(2.6691061) q[3];
sx q[3];
rz(-1.5269591) q[3];
sx q[3];
rz(2.2842177) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.040314019) q[0];
sx q[0];
rz(-0.89656985) q[0];
sx q[0];
rz(-2.8033946) q[0];
rz(-1.2942554) q[1];
sx q[1];
rz(-1.5816403) q[1];
sx q[1];
rz(-2.6370874) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0230334) q[0];
sx q[0];
rz(-1.0414062) q[0];
sx q[0];
rz(-0.27079196) q[0];
rz(0.46576969) q[2];
sx q[2];
rz(-1.2301187) q[2];
sx q[2];
rz(0.53616947) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(0.63197631) q[1];
sx q[1];
rz(-2.7512433) q[1];
sx q[1];
rz(-2.2401287) q[1];
rz(0.78318627) q[3];
sx q[3];
rz(-2.3110356) q[3];
sx q[3];
rz(-0.98379788) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(1.0304886) q[2];
sx q[2];
rz(-2.3031074) q[2];
sx q[2];
rz(-1.6476691) q[2];
rz(1.6882287) q[3];
sx q[3];
rz(-0.93795347) q[3];
sx q[3];
rz(1.7593613) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.21022739) q[0];
sx q[0];
rz(-2.2670822) q[0];
sx q[0];
rz(-2.0507623) q[0];
rz(-2.6018654) q[1];
sx q[1];
rz(-1.153774) q[1];
sx q[1];
rz(0.18879034) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.990373) q[0];
sx q[0];
rz(-1.9301206) q[0];
sx q[0];
rz(1.8381401) q[0];
rz(-pi) q[1];
rz(-1.2115057) q[2];
sx q[2];
rz(-0.28887666) q[2];
sx q[2];
rz(1.8747683) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-2.530045) q[1];
sx q[1];
rz(-1.6366742) q[1];
sx q[1];
rz(2.8987315) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.64340274) q[3];
sx q[3];
rz(-2.8426369) q[3];
sx q[3];
rz(-1.9770196) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.171689) q[2];
sx q[2];
rz(-1.0287372) q[2];
sx q[2];
rz(-1.8264654) q[2];
rz(-1.5054437) q[3];
sx q[3];
rz(-2.7841778) q[3];
sx q[3];
rz(1.2833387) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.090102) q[0];
sx q[0];
rz(-0.83054709) q[0];
sx q[0];
rz(2.8175957) q[0];
rz(-1.8404768) q[1];
sx q[1];
rz(-0.83530656) q[1];
sx q[1];
rz(-1.7623998) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1349072) q[0];
sx q[0];
rz(-2.8838257) q[0];
sx q[0];
rz(-1.2447312) q[0];
rz(-pi) q[1];
x q[1];
rz(0.46160134) q[2];
sx q[2];
rz(-1.201655) q[2];
sx q[2];
rz(2.387407) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-2.0849689) q[1];
sx q[1];
rz(-1.2722204) q[1];
sx q[1];
rz(1.741239) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.2077683) q[3];
sx q[3];
rz(-2.056042) q[3];
sx q[3];
rz(0.77321929) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(0.032701187) q[2];
sx q[2];
rz(-1.6872493) q[2];
sx q[2];
rz(2.1984055) q[2];
rz(2.8105248) q[3];
sx q[3];
rz(-1.7739242) q[3];
sx q[3];
rz(-0.29512063) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.11809764) q[0];
sx q[0];
rz(-2.5546615) q[0];
sx q[0];
rz(2.3773637) q[0];
rz(-0.14097342) q[1];
sx q[1];
rz(-0.39607513) q[1];
sx q[1];
rz(-1.2932628) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.18023597) q[0];
sx q[0];
rz(-1.1399674) q[0];
sx q[0];
rz(0.21813099) q[0];
rz(-pi) q[1];
x q[1];
rz(0.13883491) q[2];
sx q[2];
rz(-2.5359557) q[2];
sx q[2];
rz(-2.5790737) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-1.5871208) q[1];
sx q[1];
rz(-1.1677824) q[1];
sx q[1];
rz(-2.9906669) q[1];
rz(-pi) q[2];
rz(2.6929018) q[3];
sx q[3];
rz(-2.7033227) q[3];
sx q[3];
rz(-1.1718307) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(0.76688898) q[2];
sx q[2];
rz(-1.039144) q[2];
sx q[2];
rz(-0.58132201) q[2];
rz(0.86822048) q[3];
sx q[3];
rz(-0.41728443) q[3];
sx q[3];
rz(-1.8301615) q[3];
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
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.593489) q[0];
sx q[0];
rz(-2.5618401) q[0];
sx q[0];
rz(-1.2506437) q[0];
rz(-0.99682322) q[1];
sx q[1];
rz(-0.88368982) q[1];
sx q[1];
rz(-1.2876127) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8563961) q[0];
sx q[0];
rz(-2.6515793) q[0];
sx q[0];
rz(-1.6243837) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.655517) q[2];
sx q[2];
rz(-1.4176148) q[2];
sx q[2];
rz(-1.7024405) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-2.3661256) q[1];
sx q[1];
rz(-2.2062416) q[1];
sx q[1];
rz(2.8833564) q[1];
rz(-pi) q[2];
rz(3.0887293) q[3];
sx q[3];
rz(-0.82517805) q[3];
sx q[3];
rz(-2.4406274) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(1.8210956) q[2];
sx q[2];
rz(-2.7536776) q[2];
sx q[2];
rz(1.8851177) q[2];
rz(0.16658941) q[3];
sx q[3];
rz(-1.5766671) q[3];
sx q[3];
rz(-1.0725718) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2607516) q[0];
sx q[0];
rz(-2.3374225) q[0];
sx q[0];
rz(-2.8163731) q[0];
rz(-2.0064158) q[1];
sx q[1];
rz(-2.4644641) q[1];
sx q[1];
rz(0.36718711) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8242278) q[0];
sx q[0];
rz(-1.5396848) q[0];
sx q[0];
rz(-0.048006417) q[0];
rz(-pi) q[1];
x q[1];
rz(0.37819241) q[2];
sx q[2];
rz(-0.66236712) q[2];
sx q[2];
rz(0.80519245) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.90875188) q[1];
sx q[1];
rz(-0.68896657) q[1];
sx q[1];
rz(2.8244551) q[1];
x q[2];
rz(2.1807947) q[3];
sx q[3];
rz(-2.4028117) q[3];
sx q[3];
rz(1.184954) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.24511589) q[2];
sx q[2];
rz(-0.8090691) q[2];
sx q[2];
rz(2.0754576) q[2];
rz(0.079244763) q[3];
sx q[3];
rz(-0.60278046) q[3];
sx q[3];
rz(1.4762896) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
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
rz(-0.29522482) q[0];
sx q[0];
rz(-1.5548779) q[0];
sx q[0];
rz(1.5665733) q[0];
rz(-0.29905839) q[1];
sx q[1];
rz(-0.60332861) q[1];
sx q[1];
rz(0.52287846) q[1];
rz(0.085061442) q[2];
sx q[2];
rz(-0.75123514) q[2];
sx q[2];
rz(-3.0820465) q[2];
rz(2.2926033) q[3];
sx q[3];
rz(-1.6860387) q[3];
sx q[3];
rz(0.54308346) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
