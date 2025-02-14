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
rz(2.430001) q[0];
sx q[0];
rz(-1.1049668) q[0];
sx q[0];
rz(2.0018863) q[0];
rz(2.3501514) q[1];
sx q[1];
rz(-2.1005519) q[1];
sx q[1];
rz(1.1295553) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0985467) q[0];
sx q[0];
rz(-2.401863) q[0];
sx q[0];
rz(-1.3003527) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.5426251) q[2];
sx q[2];
rz(-2.2337388) q[2];
sx q[2];
rz(2.4848795) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-1.754555) q[1];
sx q[1];
rz(-0.7861775) q[1];
sx q[1];
rz(-1.992186) q[1];
rz(-pi) q[2];
rz(0.10068746) q[3];
sx q[3];
rz(-2.8071981) q[3];
sx q[3];
rz(1.7582126) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.1019885) q[2];
sx q[2];
rz(-0.81059376) q[2];
sx q[2];
rz(0.90768901) q[2];
rz(0.11861079) q[3];
sx q[3];
rz(-1.6018931) q[3];
sx q[3];
rz(2.8906726) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.94659656) q[0];
sx q[0];
rz(-0.66221607) q[0];
sx q[0];
rz(3.0372341) q[0];
rz(-1.1907499) q[1];
sx q[1];
rz(-1.2079116) q[1];
sx q[1];
rz(0.45713919) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3973244) q[0];
sx q[0];
rz(-1.5400104) q[0];
sx q[0];
rz(-1.7643614) q[0];
rz(-pi) q[1];
rz(2.2486516) q[2];
sx q[2];
rz(-1.9364076) q[2];
sx q[2];
rz(1.9302238) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-1.4310265) q[1];
sx q[1];
rz(-1.885864) q[1];
sx q[1];
rz(-0.099771413) q[1];
x q[2];
rz(0.13808226) q[3];
sx q[3];
rz(-2.2347817) q[3];
sx q[3];
rz(-0.21747227) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-2.5856058) q[2];
sx q[2];
rz(-1.1588691) q[2];
sx q[2];
rz(0.15698329) q[2];
rz(-0.5213151) q[3];
sx q[3];
rz(-2.3020404) q[3];
sx q[3];
rz(-3.1275911) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
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
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.74612015) q[0];
sx q[0];
rz(-2.5465901) q[0];
sx q[0];
rz(-2.795862) q[0];
rz(1.6861457) q[1];
sx q[1];
rz(-1.8691749) q[1];
sx q[1];
rz(1.5706496) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.81996213) q[0];
sx q[0];
rz(-1.4819711) q[0];
sx q[0];
rz(-2.1169244) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.0035946058) q[2];
sx q[2];
rz(-0.94812859) q[2];
sx q[2];
rz(-2.013157) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-1.3375255) q[1];
sx q[1];
rz(-2.9225272) q[1];
sx q[1];
rz(-2.9143224) q[1];
rz(-pi) q[2];
rz(-1.3064874) q[3];
sx q[3];
rz(-0.66625957) q[3];
sx q[3];
rz(2.8961757) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(1.7741144) q[2];
sx q[2];
rz(-0.228129) q[2];
sx q[2];
rz(0.95743123) q[2];
rz(2.0188324) q[3];
sx q[3];
rz(-1.9072396) q[3];
sx q[3];
rz(-1.3721589) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6289309) q[0];
sx q[0];
rz(-1.4224195) q[0];
sx q[0];
rz(1.5976394) q[0];
rz(1.3900025) q[1];
sx q[1];
rz(-1.3163687) q[1];
sx q[1];
rz(-1.4768627) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8771956) q[0];
sx q[0];
rz(-2.0937284) q[0];
sx q[0];
rz(-1.8967129) q[0];
rz(-pi) q[1];
rz(-1.2373015) q[2];
sx q[2];
rz(-0.25165227) q[2];
sx q[2];
rz(-2.481395) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-1.0511198) q[1];
sx q[1];
rz(-0.63544151) q[1];
sx q[1];
rz(-2.1443771) q[1];
rz(0.69227211) q[3];
sx q[3];
rz(-0.41818968) q[3];
sx q[3];
rz(1.632477) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(2.5433898) q[2];
sx q[2];
rz(-1.8566088) q[2];
sx q[2];
rz(-0.99679917) q[2];
rz(-1.8761084) q[3];
sx q[3];
rz(-2.4235453) q[3];
sx q[3];
rz(2.8640981) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0604414) q[0];
sx q[0];
rz(-1.5716946) q[0];
sx q[0];
rz(1.0720217) q[0];
rz(0.95442665) q[1];
sx q[1];
rz(-1.1773033) q[1];
sx q[1];
rz(-1.6486453) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.98534855) q[0];
sx q[0];
rz(-1.7667422) q[0];
sx q[0];
rz(-0.77134706) q[0];
rz(-1.1283595) q[2];
sx q[2];
rz(-2.2432598) q[2];
sx q[2];
rz(-2.6983124) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(1.6750329) q[1];
sx q[1];
rz(-1.7132812) q[1];
sx q[1];
rz(-0.33054466) q[1];
rz(0.29376438) q[3];
sx q[3];
rz(-2.8076594) q[3];
sx q[3];
rz(1.6373843) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-2.2299049) q[2];
sx q[2];
rz(-2.7547084) q[2];
sx q[2];
rz(1.1502728) q[2];
rz(-0.041042717) q[3];
sx q[3];
rz(-1.5951472) q[3];
sx q[3];
rz(-2.7400147) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
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
rz(0.3732442) q[0];
sx q[0];
rz(-1.1078438) q[0];
sx q[0];
rz(1.9819697) q[0];
rz(-1.4609963) q[1];
sx q[1];
rz(-0.20080876) q[1];
sx q[1];
rz(1.9083091) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4077743) q[0];
sx q[0];
rz(-0.81627699) q[0];
sx q[0];
rz(0.4305779) q[0];
x q[1];
rz(-3.0843098) q[2];
sx q[2];
rz(-1.6241239) q[2];
sx q[2];
rz(1.48207) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(1.4602051) q[1];
sx q[1];
rz(-0.80599313) q[1];
sx q[1];
rz(-0.56350033) q[1];
rz(-pi) q[2];
rz(-2.7090577) q[3];
sx q[3];
rz(-0.56946856) q[3];
sx q[3];
rz(1.7908975) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-0.50501862) q[2];
sx q[2];
rz(-0.36618149) q[2];
sx q[2];
rz(-1.9912857) q[2];
rz(-0.73927528) q[3];
sx q[3];
rz(-1.247568) q[3];
sx q[3];
rz(2.6967743) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.66733852) q[0];
sx q[0];
rz(-2.4485454) q[0];
sx q[0];
rz(2.6928103) q[0];
rz(-1.601864) q[1];
sx q[1];
rz(-1.1069143) q[1];
sx q[1];
rz(1.9805699) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0036573) q[0];
sx q[0];
rz(-1.9788392) q[0];
sx q[0];
rz(1.8613226) q[0];
x q[1];
rz(0.49893219) q[2];
sx q[2];
rz(-2.009444) q[2];
sx q[2];
rz(0.89887757) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-1.6900857) q[1];
sx q[1];
rz(-1.8645608) q[1];
sx q[1];
rz(2.7965332) q[1];
rz(-pi) q[2];
x q[2];
rz(1.9830942) q[3];
sx q[3];
rz(-2.5479043) q[3];
sx q[3];
rz(0.75324654) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-1.1686958) q[2];
sx q[2];
rz(-1.3316414) q[2];
sx q[2];
rz(-1.1807582) q[2];
rz(1.0817889) q[3];
sx q[3];
rz(-1.6094145) q[3];
sx q[3];
rz(-0.42657524) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.54843724) q[0];
sx q[0];
rz(-1.8226382) q[0];
sx q[0];
rz(-2.4539808) q[0];
rz(-1.1697191) q[1];
sx q[1];
rz(-2.1673188) q[1];
sx q[1];
rz(-2.7489472) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1448185) q[0];
sx q[0];
rz(-2.6430565) q[0];
sx q[0];
rz(1.2720435) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.4314501) q[2];
sx q[2];
rz(-1.2616488) q[2];
sx q[2];
rz(1.6935371) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(1.709014) q[1];
sx q[1];
rz(-2.6582632) q[1];
sx q[1];
rz(2.9523115) q[1];
rz(0.4057377) q[3];
sx q[3];
rz(-2.9364481) q[3];
sx q[3];
rz(2.1577378) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(0.16200599) q[2];
sx q[2];
rz(-1.0065099) q[2];
sx q[2];
rz(-1.4957734) q[2];
rz(-1.6188072) q[3];
sx q[3];
rz(-2.3706172) q[3];
sx q[3];
rz(0.60105598) q[3];
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
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.2163806) q[0];
sx q[0];
rz(-0.38327152) q[0];
sx q[0];
rz(1.8076757) q[0];
rz(1.9981617) q[1];
sx q[1];
rz(-2.3107078) q[1];
sx q[1];
rz(0.7935895) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5549907) q[0];
sx q[0];
rz(-1.7287917) q[0];
sx q[0];
rz(1.642307) q[0];
rz(-0.066448224) q[2];
sx q[2];
rz(-0.96198231) q[2];
sx q[2];
rz(-0.43283909) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.90364218) q[1];
sx q[1];
rz(-1.7667084) q[1];
sx q[1];
rz(0.51301051) q[1];
rz(-pi) q[2];
rz(-1.1203483) q[3];
sx q[3];
rz(-1.4739047) q[3];
sx q[3];
rz(-3.0894699) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(2.6952343) q[2];
sx q[2];
rz(-1.0506722) q[2];
sx q[2];
rz(-1.0515155) q[2];
rz(2.4255883) q[3];
sx q[3];
rz(-2.1456238) q[3];
sx q[3];
rz(-0.21014617) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
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
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.42651549) q[0];
sx q[0];
rz(-0.85549131) q[0];
sx q[0];
rz(2.9606384) q[0];
rz(1.1894233) q[1];
sx q[1];
rz(-1.9633429) q[1];
sx q[1];
rz(-0.98035556) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.086209379) q[0];
sx q[0];
rz(-2.6329941) q[0];
sx q[0];
rz(2.9706756) q[0];
x q[1];
rz(2.3848214) q[2];
sx q[2];
rz(-0.59322651) q[2];
sx q[2];
rz(-0.41220081) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-1.0186864) q[1];
sx q[1];
rz(-2.3748906) q[1];
sx q[1];
rz(2.8793429) q[1];
rz(-pi) q[2];
rz(1.3602363) q[3];
sx q[3];
rz(-1.9580152) q[3];
sx q[3];
rz(2.5507467) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-1.1159726) q[2];
sx q[2];
rz(-0.16416922) q[2];
sx q[2];
rz(1.3244965) q[2];
rz(-0.65274158) q[3];
sx q[3];
rz(-2.1721811) q[3];
sx q[3];
rz(3.1033707) q[3];
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
x q[0];
rz(-pi/2) q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7272335) q[0];
sx q[0];
rz(-1.4215195) q[0];
sx q[0];
rz(-0.97378578) q[0];
rz(0.7242135) q[1];
sx q[1];
rz(-1.0754633) q[1];
sx q[1];
rz(0.16574688) q[1];
rz(3.0510256) q[2];
sx q[2];
rz(-1.236785) q[2];
sx q[2];
rz(-1.0192237) q[2];
rz(-2.3694585) q[3];
sx q[3];
rz(-1.6498389) q[3];
sx q[3];
rz(-2.6966358) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
