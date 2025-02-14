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
rz(1.3266069) q[0];
sx q[0];
rz(-2.9171483) q[0];
sx q[0];
rz(1.8087968) q[0];
rz(-0.83238554) q[1];
sx q[1];
rz(4.8424911) q[1];
sx q[1];
rz(10.535156) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3334167) q[0];
sx q[0];
rz(-1.7000755) q[0];
sx q[0];
rz(2.9769633) q[0];
rz(-pi) q[1];
rz(2.2825463) q[2];
sx q[2];
rz(-1.3331604) q[2];
sx q[2];
rz(-0.29881921) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-1.8800115) q[1];
sx q[1];
rz(-3.0814897) q[1];
sx q[1];
rz(-0.48659916) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.63817231) q[3];
sx q[3];
rz(-1.320457) q[3];
sx q[3];
rz(-1.4270212) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(1.9139468) q[2];
sx q[2];
rz(-0.0089184428) q[2];
sx q[2];
rz(1.439636) q[2];
rz(-1.7281744) q[3];
sx q[3];
rz(-0.011912502) q[3];
sx q[3];
rz(0.2695151) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.444376) q[0];
sx q[0];
rz(-1.5378636) q[0];
sx q[0];
rz(-2.5268396) q[0];
rz(0.53601021) q[1];
sx q[1];
rz(-3.1159846) q[1];
sx q[1];
rz(0.33686179) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1/(6*pi)) q[0];
sx q[0];
rz(-1.4362925) q[0];
sx q[0];
rz(0.10945871) q[0];
rz(-pi) q[1];
rz(0.4763255) q[2];
sx q[2];
rz(-1.5564972) q[2];
sx q[2];
rz(-0.71029019) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-1.3276827) q[1];
sx q[1];
rz(-1.5550139) q[1];
sx q[1];
rz(-1.5266839) q[1];
x q[2];
rz(-2.9514489) q[3];
sx q[3];
rz(-1.8845673) q[3];
sx q[3];
rz(-0.47167512) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-2.2577995) q[2];
sx q[2];
rz(-0.012848583) q[2];
sx q[2];
rz(-2.4419355) q[2];
rz(0.16957016) q[3];
sx q[3];
rz(-0.01378672) q[3];
sx q[3];
rz(-0.56210303) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.70497847) q[0];
sx q[0];
rz(-0.545937) q[0];
sx q[0];
rz(2.8563232) q[0];
rz(-0.26372313) q[1];
sx q[1];
rz(-0.00058760651) q[1];
sx q[1];
rz(0.79424167) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.061210074) q[0];
sx q[0];
rz(-1.4540577) q[0];
sx q[0];
rz(-0.6878797) q[0];
rz(-pi) q[1];
x q[1];
rz(0.10767605) q[2];
sx q[2];
rz(-1.2491944) q[2];
sx q[2];
rz(1.1664101) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-1.8869739) q[1];
sx q[1];
rz(-1.5363664) q[1];
sx q[1];
rz(1.5643584) q[1];
rz(0.66624347) q[3];
sx q[3];
rz(-2.1123721) q[3];
sx q[3];
rz(0.2592087) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(0.11702015) q[2];
sx q[2];
rz(-0.046155013) q[2];
sx q[2];
rz(0.864492) q[2];
rz(-0.74508673) q[3];
sx q[3];
rz(-2.2741208) q[3];
sx q[3];
rz(0.043070506) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2660148) q[0];
sx q[0];
rz(-0.046253007) q[0];
sx q[0];
rz(0.86638802) q[0];
rz(-2.5500747) q[1];
sx q[1];
rz(-2.1951127) q[1];
sx q[1];
rz(1.0513069) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.49120228) q[0];
sx q[0];
rz(-1.5687546) q[0];
sx q[0];
rz(1.570684) q[0];
rz(1.3760304) q[2];
sx q[2];
rz(-3.1392619) q[2];
sx q[2];
rz(-2.9515206) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(1.2308553) q[1];
sx q[1];
rz(-1.9702818) q[1];
sx q[1];
rz(-0.81137701) q[1];
rz(0.71612181) q[3];
sx q[3];
rz(-2.4181626) q[3];
sx q[3];
rz(-1.0264068) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(1.5279348) q[2];
sx q[2];
rz(-3.077007) q[2];
sx q[2];
rz(-0.85395542) q[2];
rz(-2.4438786) q[3];
sx q[3];
rz(-1.8055975) q[3];
sx q[3];
rz(-2.8310988) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7455604) q[0];
sx q[0];
rz(-0.10265352) q[0];
sx q[0];
rz(2.7295617) q[0];
rz(1.8585867) q[1];
sx q[1];
rz(-0.7737776) q[1];
sx q[1];
rz(2.3903019) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.368741) q[0];
sx q[0];
rz(-1.4533459) q[0];
sx q[0];
rz(2.9020578) q[0];
rz(1.8502949) q[2];
sx q[2];
rz(-2.44254) q[2];
sx q[2];
rz(1.1115526) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-1.3739864) q[1];
sx q[1];
rz(-0.81175488) q[1];
sx q[1];
rz(2.4381258) q[1];
rz(-1.4268239) q[3];
sx q[3];
rz(-1.6917042) q[3];
sx q[3];
rz(2.6536453) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-0.46386197) q[2];
sx q[2];
rz(-0.025721392) q[2];
sx q[2];
rz(0.99979293) q[2];
rz(-2.540588) q[3];
sx q[3];
rz(-0.060636245) q[3];
sx q[3];
rz(0.84990466) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.26337013) q[0];
sx q[0];
rz(-3.0510986) q[0];
sx q[0];
rz(2.7409842) q[0];
rz(-1.406631) q[1];
sx q[1];
rz(-2.7901283) q[1];
sx q[1];
rz(-1.7469143) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.050848518) q[0];
sx q[0];
rz(-1.5904477) q[0];
sx q[0];
rz(-3.1369563) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.3371295) q[2];
sx q[2];
rz(-1.7906252) q[2];
sx q[2];
rz(1.9958082) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.7813599) q[1];
sx q[1];
rz(-2.951024) q[1];
sx q[1];
rz(3.1182454) q[1];
x q[2];
rz(-2.4056099) q[3];
sx q[3];
rz(-1.6998359) q[3];
sx q[3];
rz(2.1202212) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-2.3957735) q[2];
sx q[2];
rz(-0.58818156) q[2];
sx q[2];
rz(-2.0664717) q[2];
rz(2.6221258) q[3];
sx q[3];
rz(-2.963701) q[3];
sx q[3];
rz(1.547267) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9290685) q[0];
sx q[0];
rz(-1.2146177) q[0];
sx q[0];
rz(-1.529083) q[0];
rz(0.56135881) q[1];
sx q[1];
rz(-0.012265597) q[1];
sx q[1];
rz(-2.5711109) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9425526) q[0];
sx q[0];
rz(-1.4496452) q[0];
sx q[0];
rz(3.0899307) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.95593234) q[2];
sx q[2];
rz(-2.758965) q[2];
sx q[2];
rz(2.0418233) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(1.39064) q[1];
sx q[1];
rz(-1.5696073) q[1];
sx q[1];
rz(1.5552862) q[1];
rz(0.53855207) q[3];
sx q[3];
rz(-2.2411514) q[3];
sx q[3];
rz(0.30645257) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-3.0824215) q[2];
sx q[2];
rz(-0.26796451) q[2];
sx q[2];
rz(-1.1632261) q[2];
rz(-1.5393114) q[3];
sx q[3];
rz(-0.10207615) q[3];
sx q[3];
rz(0.99360895) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7546643) q[0];
sx q[0];
rz(-0.29351497) q[0];
sx q[0];
rz(-0.64025229) q[0];
rz(2.2786268) q[1];
sx q[1];
rz(-0.14185618) q[1];
sx q[1];
rz(2.364025) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8391188) q[0];
sx q[0];
rz(-2.5632052) q[0];
sx q[0];
rz(-0.79721398) q[0];
rz(-2.8993494) q[2];
sx q[2];
rz(-0.7787593) q[2];
sx q[2];
rz(-2.3295516) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-1.150084) q[1];
sx q[1];
rz(-1.5593632) q[1];
sx q[1];
rz(-1.4978133) q[1];
rz(2.0361774) q[3];
sx q[3];
rz(-2.0267365) q[3];
sx q[3];
rz(2.3400644) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(2.5295279) q[2];
sx q[2];
rz(-0.42403388) q[2];
sx q[2];
rz(2.6112134) q[2];
rz(-0.027675962) q[3];
sx q[3];
rz(-3.0961302) q[3];
sx q[3];
rz(1.8224705) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.46488047) q[0];
sx q[0];
rz(-0.16093971) q[0];
sx q[0];
rz(-0.93223923) q[0];
rz(-1.5981916) q[1];
sx q[1];
rz(-0.80359572) q[1];
sx q[1];
rz(-0.34206051) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.31158221) q[0];
sx q[0];
rz(-1.5556635) q[0];
sx q[0];
rz(-0.084447817) q[0];
rz(-pi) q[1];
x q[1];
rz(1.507296) q[2];
sx q[2];
rz(-0.8339774) q[2];
sx q[2];
rz(0.33499107) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-1.2310357) q[1];
sx q[1];
rz(-1.0536095) q[1];
sx q[1];
rz(-1.3523577) q[1];
rz(-pi) q[2];
rz(-0.016870528) q[3];
sx q[3];
rz(-1.3417435) q[3];
sx q[3];
rz(-3.0927174) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(3.0590234) q[2];
sx q[2];
rz(-3.1408568) q[2];
sx q[2];
rz(1.9334582) q[2];
rz(-0.9515323) q[3];
sx q[3];
rz(-3.1335242) q[3];
sx q[3];
rz(2.4212196) q[3];
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
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.35219881) q[0];
sx q[0];
rz(-0.77704) q[0];
sx q[0];
rz(0.081789516) q[0];
rz(0.31546047) q[1];
sx q[1];
rz(-0.052736484) q[1];
sx q[1];
rz(-1.9116521) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8982809) q[0];
sx q[0];
rz(-1.4139685) q[0];
sx q[0];
rz(2.9812198) q[0];
rz(-pi) q[1];
x q[1];
rz(2.4794078) q[2];
sx q[2];
rz(-0.22394584) q[2];
sx q[2];
rz(-2.0020773) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(1.1217864) q[1];
sx q[1];
rz(-0.14827327) q[1];
sx q[1];
rz(-2.2231977) q[1];
rz(-pi) q[2];
rz(2.1020402) q[3];
sx q[3];
rz(-1.4508411) q[3];
sx q[3];
rz(0.95748544) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(0.54084593) q[2];
sx q[2];
rz(-0.019286152) q[2];
sx q[2];
rz(-1.1136327) q[2];
rz(-0.039637808) q[3];
sx q[3];
rz(-3.1318635) q[3];
sx q[3];
rz(0.59687328) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
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
rz(pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6780846) q[0];
sx q[0];
rz(-1.9263374) q[0];
sx q[0];
rz(1.3081464) q[0];
rz(0.61617638) q[1];
sx q[1];
rz(-1.0721075) q[1];
sx q[1];
rz(0.20872605) q[1];
rz(-2.2439416) q[2];
sx q[2];
rz(-2.8406526) q[2];
sx q[2];
rz(0.1546897) q[2];
rz(1.3499089) q[3];
sx q[3];
rz(-1.532423) q[3];
sx q[3];
rz(1.5929089) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
