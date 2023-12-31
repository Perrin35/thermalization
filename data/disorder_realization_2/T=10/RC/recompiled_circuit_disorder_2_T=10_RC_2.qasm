OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-1.6931273) q[0];
sx q[0];
rz(-0.52283302) q[0];
sx q[0];
rz(-0.62358207) q[0];
rz(0.29016718) q[1];
sx q[1];
rz(-2.4224412) q[1];
sx q[1];
rz(-0.50049385) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1815163) q[0];
sx q[0];
rz(-1.5746563) q[0];
sx q[0];
rz(0.5998248) q[0];
x q[1];
rz(-0.093437336) q[2];
sx q[2];
rz(-1.7200025) q[2];
sx q[2];
rz(2.0051533) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(1.1653633) q[1];
sx q[1];
rz(-0.86127087) q[1];
sx q[1];
rz(2.707259) q[1];
x q[2];
rz(1.8960564) q[3];
sx q[3];
rz(-0.45913011) q[3];
sx q[3];
rz(-1.9838651) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.7913251) q[2];
sx q[2];
rz(-1.9359549) q[2];
sx q[2];
rz(1.9187437) q[2];
rz(1.6932999) q[3];
sx q[3];
rz(-0.99213123) q[3];
sx q[3];
rz(-2.1712415) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7704849) q[0];
sx q[0];
rz(-1.6587695) q[0];
sx q[0];
rz(-1.0789385) q[0];
rz(-1.7547912) q[1];
sx q[1];
rz(-2.3290122) q[1];
sx q[1];
rz(-0.66545495) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.21214813) q[0];
sx q[0];
rz(-1.11709) q[0];
sx q[0];
rz(-1.2544022) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.7416679) q[2];
sx q[2];
rz(-1.8881646) q[2];
sx q[2];
rz(-0.7764118) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-1.3033777) q[1];
sx q[1];
rz(-1.8068552) q[1];
sx q[1];
rz(-2.6462376) q[1];
rz(-0.034025107) q[3];
sx q[3];
rz(-0.861654) q[3];
sx q[3];
rz(0.89108407) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(2.5370496) q[2];
sx q[2];
rz(-1.4031354) q[2];
sx q[2];
rz(0.075142168) q[2];
rz(1.6710619) q[3];
sx q[3];
rz(-1.1276779) q[3];
sx q[3];
rz(2.4501734) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.55494088) q[0];
sx q[0];
rz(-1.8692769) q[0];
sx q[0];
rz(2.3828322) q[0];
rz(1.8485908) q[1];
sx q[1];
rz(-1.2932152) q[1];
sx q[1];
rz(-1.741515) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5132719) q[0];
sx q[0];
rz(-1.3043881) q[0];
sx q[0];
rz(-0.4011641) q[0];
x q[1];
rz(-1.7064352) q[2];
sx q[2];
rz(-1.3856158) q[2];
sx q[2];
rz(0.26091012) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(0.53576614) q[1];
sx q[1];
rz(-1.8927791) q[1];
sx q[1];
rz(0.95226007) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.8758043) q[3];
sx q[3];
rz(-1.7294356) q[3];
sx q[3];
rz(0.97914417) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-2.489958) q[2];
sx q[2];
rz(-0.48214665) q[2];
sx q[2];
rz(-0.65762323) q[2];
rz(-1.970132) q[3];
sx q[3];
rz(-1.6572584) q[3];
sx q[3];
rz(1.4191779) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.79384971) q[0];
sx q[0];
rz(-2.1934953) q[0];
sx q[0];
rz(-1.6963652) q[0];
rz(1.4472648) q[1];
sx q[1];
rz(-1.4935962) q[1];
sx q[1];
rz(2.7935374) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.30663438) q[0];
sx q[0];
rz(-1.277703) q[0];
sx q[0];
rz(0.65960633) q[0];
rz(1.939417) q[2];
sx q[2];
rz(-2.7432132) q[2];
sx q[2];
rz(2.5620661) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.066612331) q[1];
sx q[1];
rz(-2.1790824) q[1];
sx q[1];
rz(-0.10886701) q[1];
rz(1.4361037) q[3];
sx q[3];
rz(-1.4105984) q[3];
sx q[3];
rz(2.687541) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.41670123) q[2];
sx q[2];
rz(-1.3092224) q[2];
sx q[2];
rz(0.42281881) q[2];
rz(0.73741284) q[3];
sx q[3];
rz(-2.335572) q[3];
sx q[3];
rz(-0.084658682) q[3];
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
rz(pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2622862) q[0];
sx q[0];
rz(-1.8286185) q[0];
sx q[0];
rz(2.6111531) q[0];
rz(-0.92492217) q[1];
sx q[1];
rz(-1.5735807) q[1];
sx q[1];
rz(-1.8431429) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2753678) q[0];
sx q[0];
rz(-0.21731649) q[0];
sx q[0];
rz(2.2989681) q[0];
x q[1];
rz(-0.60646306) q[2];
sx q[2];
rz(-1.3242553) q[2];
sx q[2];
rz(-0.75070565) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(1.4106771) q[1];
sx q[1];
rz(-1.2123322) q[1];
sx q[1];
rz(-1.6770384) q[1];
rz(-pi) q[2];
rz(2.1634444) q[3];
sx q[3];
rz(-1.1453298) q[3];
sx q[3];
rz(-0.076171906) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.2003145) q[2];
sx q[2];
rz(-1.1555187) q[2];
sx q[2];
rz(0.47362622) q[2];
rz(-3.04223) q[3];
sx q[3];
rz(-1.2800346) q[3];
sx q[3];
rz(2.30106) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.50399238) q[0];
sx q[0];
rz(-1.3562599) q[0];
sx q[0];
rz(0.9978869) q[0];
rz(-0.87431327) q[1];
sx q[1];
rz(-1.0214146) q[1];
sx q[1];
rz(0.46674892) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3335553) q[0];
sx q[0];
rz(-1.3555962) q[0];
sx q[0];
rz(-2.4654885) q[0];
rz(-pi) q[1];
rz(1.3495965) q[2];
sx q[2];
rz(-1.1525407) q[2];
sx q[2];
rz(-2.7115371) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(1.4371498) q[1];
sx q[1];
rz(-1.4198507) q[1];
sx q[1];
rz(-1.5228065) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.7344597) q[3];
sx q[3];
rz(-1.4758849) q[3];
sx q[3];
rz(1.6223736) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.85990396) q[2];
sx q[2];
rz(-2.6624661) q[2];
sx q[2];
rz(1.5647282) q[2];
rz(-0.62670296) q[3];
sx q[3];
rz(-1.3650711) q[3];
sx q[3];
rz(2.4600162) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3307813) q[0];
sx q[0];
rz(-2.2450876) q[0];
sx q[0];
rz(-2.7217641) q[0];
rz(-0.22142521) q[1];
sx q[1];
rz(-2.6629993) q[1];
sx q[1];
rz(1.5484757) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4683068) q[0];
sx q[0];
rz(-1.5791897) q[0];
sx q[0];
rz(-0.087455672) q[0];
x q[1];
rz(-1.8553472) q[2];
sx q[2];
rz(-1.7562521) q[2];
sx q[2];
rz(0.68765771) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-0.49861136) q[1];
sx q[1];
rz(-1.6763408) q[1];
sx q[1];
rz(-1.5555744) q[1];
x q[2];
rz(-1.1057304) q[3];
sx q[3];
rz(-1.2601868) q[3];
sx q[3];
rz(-0.066699337) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(0.55591136) q[2];
sx q[2];
rz(-0.53148091) q[2];
sx q[2];
rz(-1.4038203) q[2];
rz(0.79706556) q[3];
sx q[3];
rz(-0.46204391) q[3];
sx q[3];
rz(1.2169303) q[3];
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
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4306915) q[0];
sx q[0];
rz(-0.71390188) q[0];
sx q[0];
rz(-2.8523493) q[0];
rz(-0.62492433) q[1];
sx q[1];
rz(-1.9754675) q[1];
sx q[1];
rz(-1.3141059) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8970012) q[0];
sx q[0];
rz(-2.3750711) q[0];
sx q[0];
rz(-0.15890973) q[0];
rz(-0.74771379) q[2];
sx q[2];
rz(-1.3249825) q[2];
sx q[2];
rz(-2.4112548) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(0.18409477) q[1];
sx q[1];
rz(-0.1322608) q[1];
sx q[1];
rz(-1.9930507) q[1];
rz(2.3280011) q[3];
sx q[3];
rz(-1.3435257) q[3];
sx q[3];
rz(-3.1384625) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-0.73359314) q[2];
sx q[2];
rz(-1.5635798) q[2];
sx q[2];
rz(-1.4286263) q[2];
rz(-2.1777878) q[3];
sx q[3];
rz(-2.0644085) q[3];
sx q[3];
rz(-1.0296286) q[3];
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
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6190417) q[0];
sx q[0];
rz(-1.4551117) q[0];
sx q[0];
rz(1.8956986) q[0];
rz(-3.030581) q[1];
sx q[1];
rz(-1.9440034) q[1];
sx q[1];
rz(-0.54661173) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3205991) q[0];
sx q[0];
rz(-0.59251596) q[0];
sx q[0];
rz(-0.86235637) q[0];
rz(-pi) q[1];
x q[1];
rz(0.75066363) q[2];
sx q[2];
rz(-2.6235136) q[2];
sx q[2];
rz(0.93875611) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-2.1281801) q[1];
sx q[1];
rz(-1.9617404) q[1];
sx q[1];
rz(0.18545111) q[1];
rz(-pi) q[2];
x q[2];
rz(2.2253753) q[3];
sx q[3];
rz(-1.9700053) q[3];
sx q[3];
rz(0.44772128) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-1.0299915) q[2];
sx q[2];
rz(-2.2932055) q[2];
sx q[2];
rz(-1.4477504) q[2];
rz(2.0041806) q[3];
sx q[3];
rz(-1.3237938) q[3];
sx q[3];
rz(0.64731961) q[3];
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
sx q[0];
rz(-pi) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.85957134) q[0];
sx q[0];
rz(-2.8347926) q[0];
sx q[0];
rz(-2.4243673) q[0];
rz(1.2099129) q[1];
sx q[1];
rz(-2.8094493) q[1];
sx q[1];
rz(-0.70770121) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.74117888) q[0];
sx q[0];
rz(-0.78010633) q[0];
sx q[0];
rz(-2.0692503) q[0];
rz(-pi) q[1];
x q[1];
rz(1.0102347) q[2];
sx q[2];
rz(-2.8786504) q[2];
sx q[2];
rz(1.5514785) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-2.6245539) q[1];
sx q[1];
rz(-2.3561764) q[1];
sx q[1];
rz(-2.024827) q[1];
rz(-2.7865949) q[3];
sx q[3];
rz(-0.98155752) q[3];
sx q[3];
rz(-2.0139351) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.5132961) q[2];
sx q[2];
rz(-2.4847023) q[2];
sx q[2];
rz(0.90325242) q[2];
rz(-1.5385657) q[3];
sx q[3];
rz(-2.2730946) q[3];
sx q[3];
rz(-0.85047754) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
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
rz(-0.83508867) q[0];
sx q[0];
rz(-2.7764414) q[0];
sx q[0];
rz(2.2055702) q[0];
rz(-0.8159591) q[1];
sx q[1];
rz(-0.42146704) q[1];
sx q[1];
rz(-2.0889919) q[1];
rz(-2.6782398) q[2];
sx q[2];
rz(-1.9909161) q[2];
sx q[2];
rz(-1.9009895) q[2];
rz(-1.5296616) q[3];
sx q[3];
rz(-1.516468) q[3];
sx q[3];
rz(-0.80249912) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
