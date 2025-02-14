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
rz(2.6142081) q[0];
sx q[0];
rz(-2.5684147) q[0];
sx q[0];
rz(2.5249124) q[0];
rz(-1.0083899) q[1];
sx q[1];
rz(3.8808793) q[1];
sx q[1];
rz(9.0864658) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.2074006) q[0];
sx q[0];
rz(-1.3284042) q[0];
sx q[0];
rz(1.6079748) q[0];
x q[1];
rz(0.7912999) q[2];
sx q[2];
rz(-1.2630579) q[2];
sx q[2];
rz(1.5200293) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-1.4614726) q[1];
sx q[1];
rz(-2.1109889) q[1];
sx q[1];
rz(-1.4137181) q[1];
x q[2];
rz(2.9659392) q[3];
sx q[3];
rz(-0.82161108) q[3];
sx q[3];
rz(2.780811) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(1.3931291) q[2];
sx q[2];
rz(-1.4145565) q[2];
sx q[2];
rz(0.65790042) q[2];
rz(-2.6864478) q[3];
sx q[3];
rz(-0.15076605) q[3];
sx q[3];
rz(1.8337102) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8667792) q[0];
sx q[0];
rz(-1.4115189) q[0];
sx q[0];
rz(-2.859512) q[0];
rz(0.2433978) q[1];
sx q[1];
rz(-1.2671821) q[1];
sx q[1];
rz(0.74877053) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.98784763) q[0];
sx q[0];
rz(-1.7637588) q[0];
sx q[0];
rz(0.70867507) q[0];
x q[1];
rz(-1.5087328) q[2];
sx q[2];
rz(-0.12189874) q[2];
sx q[2];
rz(-0.97165194) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(1.1324596) q[1];
sx q[1];
rz(-1.0693502) q[1];
sx q[1];
rz(-3.0363068) q[1];
rz(-pi) q[2];
rz(0.19939662) q[3];
sx q[3];
rz(-1.6391338) q[3];
sx q[3];
rz(2.0082366) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(2.8613209) q[2];
sx q[2];
rz(-1.5254285) q[2];
sx q[2];
rz(-1.8005499) q[2];
rz(1.9848112) q[3];
sx q[3];
rz(-2.3564434) q[3];
sx q[3];
rz(2.3780499) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
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
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.84083104) q[0];
sx q[0];
rz(-2.9756727) q[0];
sx q[0];
rz(2.4842343) q[0];
rz(-1.1454469) q[1];
sx q[1];
rz(-0.95445389) q[1];
sx q[1];
rz(-2.3232536) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9377146) q[0];
sx q[0];
rz(-2.7174207) q[0];
sx q[0];
rz(0.61221497) q[0];
x q[1];
rz(-3.0938593) q[2];
sx q[2];
rz(-1.2558508) q[2];
sx q[2];
rz(2.3746109) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.044825252) q[1];
sx q[1];
rz(-0.92508537) q[1];
sx q[1];
rz(0.94618465) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.0243689) q[3];
sx q[3];
rz(-1.6370341) q[3];
sx q[3];
rz(0.57285786) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(2.1546617) q[2];
sx q[2];
rz(-1.2620474) q[2];
sx q[2];
rz(-1.7884802) q[2];
rz(2.1766369) q[3];
sx q[3];
rz(-1.302364) q[3];
sx q[3];
rz(-2.7195948) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.44593909) q[0];
sx q[0];
rz(-0.66534477) q[0];
sx q[0];
rz(1.2009784) q[0];
rz(0.0032084223) q[1];
sx q[1];
rz(-1.4326347) q[1];
sx q[1];
rz(0.23689717) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.1334522) q[0];
sx q[0];
rz(-2.4137375) q[0];
sx q[0];
rz(-0.86489622) q[0];
rz(-pi) q[1];
x q[1];
rz(2.5043244) q[2];
sx q[2];
rz(-1.6760525) q[2];
sx q[2];
rz(2.1861475) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.8767406) q[1];
sx q[1];
rz(-2.5036466) q[1];
sx q[1];
rz(-1.9015584) q[1];
rz(-pi) q[2];
rz(2.3637216) q[3];
sx q[3];
rz(-2.155455) q[3];
sx q[3];
rz(-2.8608222) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(1.0995471) q[2];
sx q[2];
rz(-1.8501661) q[2];
sx q[2];
rz(2.715204) q[2];
rz(2.1060627) q[3];
sx q[3];
rz(-1.6798881) q[3];
sx q[3];
rz(0.6216014) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.34430382) q[0];
sx q[0];
rz(-3.099589) q[0];
sx q[0];
rz(2.7916743) q[0];
rz(1.2087076) q[1];
sx q[1];
rz(-1.8588926) q[1];
sx q[1];
rz(-0.0016317687) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.37536538) q[0];
sx q[0];
rz(-1.3266801) q[0];
sx q[0];
rz(1.3993511) q[0];
rz(-pi) q[1];
x q[1];
rz(2.8441378) q[2];
sx q[2];
rz(-1.607873) q[2];
sx q[2];
rz(-0.86264474) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(0.35910159) q[1];
sx q[1];
rz(-1.3939855) q[1];
sx q[1];
rz(0.048080877) q[1];
rz(-0.91301133) q[3];
sx q[3];
rz(-0.33644852) q[3];
sx q[3];
rz(-0.89321619) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(2.1899167) q[2];
sx q[2];
rz(-2.2168171) q[2];
sx q[2];
rz(2.9803989) q[2];
rz(-0.4392043) q[3];
sx q[3];
rz(-2.3930211) q[3];
sx q[3];
rz(-2.2883033) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.31140232) q[0];
sx q[0];
rz(-1.4051733) q[0];
sx q[0];
rz(-0.79469529) q[0];
rz(1.3658124) q[1];
sx q[1];
rz(-0.8005442) q[1];
sx q[1];
rz(-0.47439233) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.2524388) q[0];
sx q[0];
rz(-1.948101) q[0];
sx q[0];
rz(-2.120156) q[0];
rz(-pi) q[1];
x q[1];
rz(2.5849336) q[2];
sx q[2];
rz(-1.2534598) q[2];
sx q[2];
rz(-2.5450626) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(1.5833334) q[1];
sx q[1];
rz(-0.72860241) q[1];
sx q[1];
rz(2.5069951) q[1];
rz(-pi) q[2];
rz(-2.1544564) q[3];
sx q[3];
rz(-1.2452092) q[3];
sx q[3];
rz(-1.1914355) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(0.61364335) q[2];
sx q[2];
rz(-1.8972634) q[2];
sx q[2];
rz(-0.5274241) q[2];
rz(0.6066277) q[3];
sx q[3];
rz(-0.18692034) q[3];
sx q[3];
rz(-2.0691779) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.2622751) q[0];
sx q[0];
rz(-0.19392218) q[0];
sx q[0];
rz(0.79757565) q[0];
rz(2.2480615) q[1];
sx q[1];
rz(-0.83502665) q[1];
sx q[1];
rz(2.8614047) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2036616) q[0];
sx q[0];
rz(-1.1755921) q[0];
sx q[0];
rz(-2.9523938) q[0];
x q[1];
rz(-2.1970388) q[2];
sx q[2];
rz(-1.0172067) q[2];
sx q[2];
rz(-1.1819201) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-0.60498991) q[1];
sx q[1];
rz(-2.0532673) q[1];
sx q[1];
rz(-3.1274113) q[1];
rz(0.35983054) q[3];
sx q[3];
rz(-2.1323279) q[3];
sx q[3];
rz(2.7240745) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(1.7402652) q[2];
sx q[2];
rz(-1.3046616) q[2];
sx q[2];
rz(-1.2124445) q[2];
rz(0.23981833) q[3];
sx q[3];
rz(-0.66767728) q[3];
sx q[3];
rz(-0.56813204) q[3];
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
rz(2.3577393) q[0];
sx q[0];
rz(-1.724406) q[0];
sx q[0];
rz(1.3215815) q[0];
rz(-2.9603738) q[1];
sx q[1];
rz(-1.3316589) q[1];
sx q[1];
rz(-0.063057335) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2305166) q[0];
sx q[0];
rz(-2.3421093) q[0];
sx q[0];
rz(1.37079) q[0];
rz(-pi) q[1];
rz(-0.087788344) q[2];
sx q[2];
rz(-1.8551747) q[2];
sx q[2];
rz(1.2149917) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-1.6136421) q[1];
sx q[1];
rz(-0.41955742) q[1];
sx q[1];
rz(2.502524) q[1];
x q[2];
rz(-1.5884871) q[3];
sx q[3];
rz(-2.9790386) q[3];
sx q[3];
rz(-1.7944195) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(2.7276089) q[2];
sx q[2];
rz(-2.0435645) q[2];
sx q[2];
rz(1.46924) q[2];
rz(1.3937048) q[3];
sx q[3];
rz(-1.7017476) q[3];
sx q[3];
rz(-2.9318455) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
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
rz(-2.6032228) q[0];
sx q[0];
rz(-2.5312238) q[0];
sx q[0];
rz(0.99697733) q[0];
rz(2.3485377) q[1];
sx q[1];
rz(-1.4184364) q[1];
sx q[1];
rz(2.0319895) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.0080896688) q[0];
sx q[0];
rz(-0.64787302) q[0];
sx q[0];
rz(-2.6459207) q[0];
rz(-pi) q[1];
x q[1];
rz(2.419554) q[2];
sx q[2];
rz(-1.2424316) q[2];
sx q[2];
rz(-1.9313081) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.80191699) q[1];
sx q[1];
rz(-1.1028506) q[1];
sx q[1];
rz(0.94593327) q[1];
x q[2];
rz(1.7392735) q[3];
sx q[3];
rz(-0.94120294) q[3];
sx q[3];
rz(1.633267) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-0.41254607) q[2];
sx q[2];
rz(-1.46526) q[2];
sx q[2];
rz(1.868978) q[2];
rz(-2.0160969) q[3];
sx q[3];
rz(-2.9477305) q[3];
sx q[3];
rz(0.079843609) q[3];
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
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7703055) q[0];
sx q[0];
rz(-1.9885539) q[0];
sx q[0];
rz(-3.0433997) q[0];
rz(-1.498361) q[1];
sx q[1];
rz(-1.9941092) q[1];
sx q[1];
rz(0.8383382) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6107008) q[0];
sx q[0];
rz(-2.3641788) q[0];
sx q[0];
rz(-2.8882508) q[0];
rz(-pi) q[1];
rz(-0.21282332) q[2];
sx q[2];
rz(-2.8838919) q[2];
sx q[2];
rz(2.3824196) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(0.34217087) q[1];
sx q[1];
rz(-1.8525181) q[1];
sx q[1];
rz(-1.1200302) q[1];
rz(0.43180744) q[3];
sx q[3];
rz(-1.7718959) q[3];
sx q[3];
rz(-2.3922298) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(1.3457634) q[2];
sx q[2];
rz(-0.886262) q[2];
sx q[2];
rz(2.8161827) q[2];
rz(2.365153) q[3];
sx q[3];
rz(-1.0151981) q[3];
sx q[3];
rz(-2.5344892) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0008739) q[0];
sx q[0];
rz(-1.9372531) q[0];
sx q[0];
rz(-1.0304864) q[0];
rz(-2.8554032) q[1];
sx q[1];
rz(-1.20594) q[1];
sx q[1];
rz(1.1084569) q[1];
rz(0.76079615) q[2];
sx q[2];
rz(-1.8626871) q[2];
sx q[2];
rz(-2.5260851) q[2];
rz(0.48578942) q[3];
sx q[3];
rz(-0.64668568) q[3];
sx q[3];
rz(0.034737094) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
