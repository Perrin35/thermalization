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
rz(-2.0120373) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4576042) q[0];
sx q[0];
rz(-0.86375551) q[0];
sx q[0];
rz(0.23914214) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.036058124) q[2];
sx q[2];
rz(-2.4781422) q[2];
sx q[2];
rz(-0.61095881) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-1.1892231) q[1];
sx q[1];
rz(-0.86878759) q[1];
sx q[1];
rz(-2.7527806) q[1];
rz(-pi) q[2];
rz(-1.6057062) q[3];
sx q[3];
rz(-1.9034317) q[3];
sx q[3];
rz(-1.2768318) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(2.1019885) q[2];
sx q[2];
rz(-2.3309989) q[2];
sx q[2];
rz(-0.90768901) q[2];
rz(0.11861079) q[3];
sx q[3];
rz(-1.5396996) q[3];
sx q[3];
rz(-2.8906726) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
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
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.94659656) q[0];
sx q[0];
rz(-2.4793766) q[0];
sx q[0];
rz(-0.10435852) q[0];
rz(1.1907499) q[1];
sx q[1];
rz(-1.2079116) q[1];
sx q[1];
rz(2.6844535) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.16743827) q[0];
sx q[0];
rz(-1.3773241) q[0];
sx q[0];
rz(-0.03137145) q[0];
x q[1];
rz(2.2486516) q[2];
sx q[2];
rz(-1.9364076) q[2];
sx q[2];
rz(-1.2113689) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-3.0328331) q[1];
sx q[1];
rz(-1.4759513) q[1];
sx q[1];
rz(1.2542568) q[1];
rz(-3.0035104) q[3];
sx q[3];
rz(-0.90681091) q[3];
sx q[3];
rz(0.21747227) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-2.5856058) q[2];
sx q[2];
rz(-1.1588691) q[2];
sx q[2];
rz(0.15698329) q[2];
rz(-2.6202776) q[3];
sx q[3];
rz(-2.3020404) q[3];
sx q[3];
rz(-0.014001525) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.74612015) q[0];
sx q[0];
rz(-0.59500256) q[0];
sx q[0];
rz(2.795862) q[0];
rz(-1.6861457) q[1];
sx q[1];
rz(-1.8691749) q[1];
sx q[1];
rz(-1.5706496) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.60590505) q[0];
sx q[0];
rz(-2.5890124) q[0];
sx q[0];
rz(-1.4009757) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.0035946058) q[2];
sx q[2];
rz(-2.1934641) q[2];
sx q[2];
rz(2.013157) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(1.8040672) q[1];
sx q[1];
rz(-2.9225272) q[1];
sx q[1];
rz(0.22727025) q[1];
rz(1.3064874) q[3];
sx q[3];
rz(-0.66625957) q[3];
sx q[3];
rz(-2.8961757) q[3];
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
rz(1.1227603) q[3];
sx q[3];
rz(-1.2343531) q[3];
sx q[3];
rz(1.7694337) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.51266176) q[0];
sx q[0];
rz(-1.4224195) q[0];
sx q[0];
rz(1.5439532) q[0];
rz(-1.7515901) q[1];
sx q[1];
rz(-1.3163687) q[1];
sx q[1];
rz(-1.4768627) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.66946736) q[0];
sx q[0];
rz(-0.60807121) q[0];
sx q[0];
rz(-2.6345992) q[0];
rz(-pi) q[1];
x q[1];
rz(1.2373015) q[2];
sx q[2];
rz(-0.25165227) q[2];
sx q[2];
rz(-0.66019765) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.0904729) q[1];
sx q[1];
rz(-2.5061511) q[1];
sx q[1];
rz(2.1443771) q[1];
rz(-1.847193) q[3];
sx q[3];
rz(-1.8887465) q[3];
sx q[3];
rz(-0.89561392) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.5433898) q[2];
sx q[2];
rz(-1.8566088) q[2];
sx q[2];
rz(2.1447935) q[2];
rz(-1.2654842) q[3];
sx q[3];
rz(-0.71804738) q[3];
sx q[3];
rz(-0.27749458) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(pi/2) q[2];
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
rz(-1.9642893) q[1];
sx q[1];
rz(-1.4929474) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.98534855) q[0];
sx q[0];
rz(-1.3748504) q[0];
sx q[0];
rz(-0.77134706) q[0];
rz(-pi) q[1];
x q[1];
rz(2.6482539) q[2];
sx q[2];
rz(-0.78560053) q[2];
sx q[2];
rz(-1.0935022) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-0.055549057) q[1];
sx q[1];
rz(-1.2437268) q[1];
sx q[1];
rz(-1.4202761) q[1];
rz(1.6709153) q[3];
sx q[3];
rz(-1.2516875) q[3];
sx q[3];
rz(-1.1943194) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-0.91168779) q[2];
sx q[2];
rz(-2.7547084) q[2];
sx q[2];
rz(1.9913199) q[2];
rz(3.1005499) q[3];
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
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7683485) q[0];
sx q[0];
rz(-1.1078438) q[0];
sx q[0];
rz(-1.1596229) q[0];
rz(-1.4609963) q[1];
sx q[1];
rz(-2.9407839) q[1];
sx q[1];
rz(-1.9083091) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4077743) q[0];
sx q[0];
rz(-0.81627699) q[0];
sx q[0];
rz(-0.4305779) q[0];
rz(-pi) q[1];
rz(1.6242113) q[2];
sx q[2];
rz(-1.513595) q[2];
sx q[2];
rz(3.055923) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(2.199904) q[1];
sx q[1];
rz(-2.2268128) q[1];
sx q[1];
rz(1.0628878) q[1];
x q[2];
rz(-2.7090577) q[3];
sx q[3];
rz(-2.5721241) q[3];
sx q[3];
rz(-1.7908975) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-0.50501862) q[2];
sx q[2];
rz(-0.36618149) q[2];
sx q[2];
rz(-1.9912857) q[2];
rz(0.73927528) q[3];
sx q[3];
rz(-1.8940247) q[3];
sx q[3];
rz(2.6967743) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
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
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.66733852) q[0];
sx q[0];
rz(-0.69304729) q[0];
sx q[0];
rz(2.6928103) q[0];
rz(-1.5397286) q[1];
sx q[1];
rz(-1.1069143) q[1];
sx q[1];
rz(-1.9805699) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5509508) q[0];
sx q[0];
rz(-1.30473) q[0];
sx q[0];
rz(-0.42386415) q[0];
x q[1];
rz(-0.77552253) q[2];
sx q[2];
rz(-2.4897414) q[2];
sx q[2];
rz(-1.8076123) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(2.9185982) q[1];
sx q[1];
rz(-1.9004993) q[1];
sx q[1];
rz(-1.2597643) q[1];
rz(-2.8774639) q[3];
sx q[3];
rz(-2.1089411) q[3];
sx q[3];
rz(1.9028185) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(1.1686958) q[2];
sx q[2];
rz(-1.3316414) q[2];
sx q[2];
rz(-1.9608344) q[2];
rz(-1.0817889) q[3];
sx q[3];
rz(-1.5321782) q[3];
sx q[3];
rz(-0.42657524) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5931554) q[0];
sx q[0];
rz(-1.8226382) q[0];
sx q[0];
rz(2.4539808) q[0];
rz(1.1697191) q[1];
sx q[1];
rz(-2.1673188) q[1];
sx q[1];
rz(-0.39264548) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3340281) q[0];
sx q[0];
rz(-1.0962209) q[0];
sx q[0];
rz(2.9827098) q[0];
x q[1];
rz(-2.7313825) q[2];
sx q[2];
rz(-2.8034119) q[2];
sx q[2];
rz(1.2615801) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-3.1117797) q[1];
sx q[1];
rz(-1.6583484) q[1];
sx q[1];
rz(0.47595166) q[1];
rz(-pi) q[2];
rz(-1.6527376) q[3];
sx q[3];
rz(-1.382516) q[3];
sx q[3];
rz(0.57037607) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.9795867) q[2];
sx q[2];
rz(-2.1350828) q[2];
sx q[2];
rz(-1.4957734) q[2];
rz(1.6188072) q[3];
sx q[3];
rz(-2.3706172) q[3];
sx q[3];
rz(-0.60105598) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.2163806) q[0];
sx q[0];
rz(-2.7583211) q[0];
sx q[0];
rz(1.8076757) q[0];
rz(1.9981617) q[1];
sx q[1];
rz(-2.3107078) q[1];
sx q[1];
rz(0.7935895) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1593678) q[0];
sx q[0];
rz(-0.17330226) q[0];
sx q[0];
rz(-2.7200219) q[0];
rz(-pi) q[1];
x q[1];
rz(2.1806472) q[2];
sx q[2];
rz(-1.5163002) q[2];
sx q[2];
rz(1.1759963) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-2.1417228) q[1];
sx q[1];
rz(-2.595584) q[1];
sx q[1];
rz(2.757339) q[1];
rz(-1.7904441) q[3];
sx q[3];
rz(-2.6815412) q[3];
sx q[3];
rz(-1.3212412) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-0.44635832) q[2];
sx q[2];
rz(-1.0506722) q[2];
sx q[2];
rz(-1.0515155) q[2];
rz(-0.71600437) q[3];
sx q[3];
rz(-0.99596888) q[3];
sx q[3];
rz(-2.9314465) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7150772) q[0];
sx q[0];
rz(-2.2861013) q[0];
sx q[0];
rz(-0.18095428) q[0];
rz(-1.1894233) q[1];
sx q[1];
rz(-1.9633429) q[1];
sx q[1];
rz(0.98035556) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.086209379) q[0];
sx q[0];
rz(-0.50859857) q[0];
sx q[0];
rz(2.9706756) q[0];
rz(0.75677121) q[2];
sx q[2];
rz(-0.59322651) q[2];
sx q[2];
rz(0.41220081) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.1229062) q[1];
sx q[1];
rz(-0.76670206) q[1];
sx q[1];
rz(-0.26224978) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.7465024) q[3];
sx q[3];
rz(-1.7655585) q[3];
sx q[3];
rz(0.89941809) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.1159726) q[2];
sx q[2];
rz(-0.16416922) q[2];
sx q[2];
rz(-1.8170961) q[2];
rz(-0.65274158) q[3];
sx q[3];
rz(-2.1721811) q[3];
sx q[3];
rz(-0.038221922) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.41435913) q[0];
sx q[0];
rz(-1.7200732) q[0];
sx q[0];
rz(2.1678069) q[0];
rz(-2.4173792) q[1];
sx q[1];
rz(-1.0754633) q[1];
sx q[1];
rz(0.16574688) q[1];
rz(-3.0510256) q[2];
sx q[2];
rz(-1.9048077) q[2];
sx q[2];
rz(2.122369) q[2];
rz(0.11304819) q[3];
sx q[3];
rz(-2.3662576) q[3];
sx q[3];
rz(1.934847) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
