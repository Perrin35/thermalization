OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(2.9053469) q[0];
sx q[0];
rz(-0.72609225) q[0];
sx q[0];
rz(-0.2015764) q[0];
rz(-2.6456614) q[1];
sx q[1];
rz(-2.6013241) q[1];
sx q[1];
rz(0.93710605) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1256589) q[0];
sx q[0];
rz(-2.4056245) q[0];
sx q[0];
rz(-2.486881) q[0];
x q[1];
rz(2.063077) q[2];
sx q[2];
rz(-1.0168076) q[2];
sx q[2];
rz(0.89821494) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-1.1569251) q[1];
sx q[1];
rz(-2.942454) q[1];
sx q[1];
rz(-1.9763293) q[1];
rz(-pi) q[2];
rz(-1.9107781) q[3];
sx q[3];
rz(-0.28453207) q[3];
sx q[3];
rz(-1.1839379) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(0.15930882) q[2];
sx q[2];
rz(-2.8757877) q[2];
sx q[2];
rz(0.086159555) q[2];
rz(-2.384095) q[3];
sx q[3];
rz(-2.3870654) q[3];
sx q[3];
rz(-1.0936201) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
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
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5646097) q[0];
sx q[0];
rz(-1.4819205) q[0];
sx q[0];
rz(1.8151059) q[0];
rz(-1.2558698) q[1];
sx q[1];
rz(-1.565226) q[1];
sx q[1];
rz(2.870141) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.87646987) q[0];
sx q[0];
rz(-1.0713097) q[0];
sx q[0];
rz(-0.45544099) q[0];
rz(-pi) q[1];
rz(-0.51606744) q[2];
sx q[2];
rz(-1.3000254) q[2];
sx q[2];
rz(-2.5004636) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(0.27122341) q[1];
sx q[1];
rz(-1.9454114) q[1];
sx q[1];
rz(-2.3323374) q[1];
x q[2];
rz(1.9430964) q[3];
sx q[3];
rz(-0.79685235) q[3];
sx q[3];
rz(0.70355319) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(2.0469971) q[2];
sx q[2];
rz(-1.7589898) q[2];
sx q[2];
rz(-0.57717741) q[2];
rz(0.92352891) q[3];
sx q[3];
rz(-2.2369592) q[3];
sx q[3];
rz(1.7318168) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
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
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.24401027) q[0];
sx q[0];
rz(-2.0799461) q[0];
sx q[0];
rz(-2.999021) q[0];
rz(1.7890731) q[1];
sx q[1];
rz(-2.1058154) q[1];
sx q[1];
rz(-0.20908633) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.81461834) q[0];
sx q[0];
rz(-2.1332392) q[0];
sx q[0];
rz(2.4787865) q[0];
rz(-1.1582583) q[2];
sx q[2];
rz(-0.87450714) q[2];
sx q[2];
rz(1.0227026) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-1.0774539) q[1];
sx q[1];
rz(-1.1541379) q[1];
sx q[1];
rz(2.0616848) q[1];
rz(-pi) q[2];
rz(2.530982) q[3];
sx q[3];
rz(-2.7770677) q[3];
sx q[3];
rz(-0.025346905) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(1.7769988) q[2];
sx q[2];
rz(-0.89407388) q[2];
sx q[2];
rz(2.7374632) q[2];
rz(-1.2858307) q[3];
sx q[3];
rz(-2.0127681) q[3];
sx q[3];
rz(-2.9742441) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[3];
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
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4362815) q[0];
sx q[0];
rz(-3.0349338) q[0];
sx q[0];
rz(0.96281111) q[0];
rz(-0.46936938) q[1];
sx q[1];
rz(-2.5517187) q[1];
sx q[1];
rz(-0.00096360047) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9199333) q[0];
sx q[0];
rz(-0.98826212) q[0];
sx q[0];
rz(0.68908738) q[0];
x q[1];
rz(-1.5702815) q[2];
sx q[2];
rz(-1.6948023) q[2];
sx q[2];
rz(-0.73881432) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-1.2767267) q[1];
sx q[1];
rz(-0.58502561) q[1];
sx q[1];
rz(1.8086955) q[1];
rz(-pi) q[2];
rz(-2.5721781) q[3];
sx q[3];
rz(-1.6462407) q[3];
sx q[3];
rz(0.30573341) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-2.0659539) q[2];
sx q[2];
rz(-1.5116189) q[2];
sx q[2];
rz(-3.0299419) q[2];
rz(-2.3305437) q[3];
sx q[3];
rz(-2.6871197) q[3];
sx q[3];
rz(3.1276935) q[3];
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
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.951293) q[0];
sx q[0];
rz(-2.2409029) q[0];
sx q[0];
rz(-0.35650373) q[0];
rz(-0.50645343) q[1];
sx q[1];
rz(-1.7849779) q[1];
sx q[1];
rz(0.26062632) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.31291744) q[0];
sx q[0];
rz(-0.14794359) q[0];
sx q[0];
rz(-0.28767985) q[0];
x q[1];
rz(2.3512958) q[2];
sx q[2];
rz(-2.7919263) q[2];
sx q[2];
rz(-0.5860354) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(1.534429) q[1];
sx q[1];
rz(-2.5669614) q[1];
sx q[1];
rz(-0.032392153) q[1];
x q[2];
rz(-0.14560933) q[3];
sx q[3];
rz(-1.7123316) q[3];
sx q[3];
rz(1.4841929) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-1.9115209) q[2];
sx q[2];
rz(-1.6610961) q[2];
sx q[2];
rz(2.9122706) q[2];
rz(0.54245943) q[3];
sx q[3];
rz(-2.8312603) q[3];
sx q[3];
rz(-2.1054161) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3689573) q[0];
sx q[0];
rz(-2.2326523) q[0];
sx q[0];
rz(1.4916346) q[0];
rz(-2.1024599) q[1];
sx q[1];
rz(-1.2960478) q[1];
sx q[1];
rz(-1.414149) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6468069) q[0];
sx q[0];
rz(-1.2959058) q[0];
sx q[0];
rz(-1.5216212) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.44712375) q[2];
sx q[2];
rz(-0.56338718) q[2];
sx q[2];
rz(-0.901957) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(0.27281877) q[1];
sx q[1];
rz(-2.0588074) q[1];
sx q[1];
rz(-2.6150319) q[1];
rz(-pi) q[2];
rz(1.9023444) q[3];
sx q[3];
rz(-2.9132531) q[3];
sx q[3];
rz(-2.6339298) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-2.002939) q[2];
sx q[2];
rz(-2.5230375) q[2];
sx q[2];
rz(-3.1138528) q[2];
rz(2.6489143) q[3];
sx q[3];
rz(-1.6511107) q[3];
sx q[3];
rz(-1.7475351) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7572927) q[0];
sx q[0];
rz(-1.5889656) q[0];
sx q[0];
rz(-3.0199155) q[0];
rz(-1.1514459) q[1];
sx q[1];
rz(-0.45184389) q[1];
sx q[1];
rz(0.40245232) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.84797317) q[0];
sx q[0];
rz(-1.8403887) q[0];
sx q[0];
rz(2.2718391) q[0];
rz(-pi) q[1];
x q[1];
rz(1.5624814) q[2];
sx q[2];
rz(-1.6168211) q[2];
sx q[2];
rz(2.7260821) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-0.35343364) q[1];
sx q[1];
rz(-2.806059) q[1];
sx q[1];
rz(-2.5787756) q[1];
rz(-pi) q[2];
rz(-1.2225371) q[3];
sx q[3];
rz(-1.5969443) q[3];
sx q[3];
rz(-2.1611283) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(3.1272614) q[2];
sx q[2];
rz(-1.1703337) q[2];
sx q[2];
rz(-0.86223117) q[2];
rz(2.6640653) q[3];
sx q[3];
rz(-1.9600441) q[3];
sx q[3];
rz(2.2235218) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.35571337) q[0];
sx q[0];
rz(-0.15545758) q[0];
sx q[0];
rz(0.73295897) q[0];
rz(-3.0015302) q[1];
sx q[1];
rz(-0.99761325) q[1];
sx q[1];
rz(-2.1070811) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2037172) q[0];
sx q[0];
rz(-2.8563742) q[0];
sx q[0];
rz(0.9271778) q[0];
rz(-1.6413692) q[2];
sx q[2];
rz(-1.4683873) q[2];
sx q[2];
rz(-2.9533353) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-2.5518783) q[1];
sx q[1];
rz(-2.1363746) q[1];
sx q[1];
rz(0.043885529) q[1];
x q[2];
rz(-0.23969527) q[3];
sx q[3];
rz(-0.42907676) q[3];
sx q[3];
rz(2.808188) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(0.90116477) q[2];
sx q[2];
rz(-1.9463836) q[2];
sx q[2];
rz(0.77587664) q[2];
rz(-2.4173229) q[3];
sx q[3];
rz(-0.80562076) q[3];
sx q[3];
rz(-1.4340713) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.6999321) q[0];
sx q[0];
rz(-2.3775546) q[0];
sx q[0];
rz(3.124776) q[0];
rz(-0.018521221) q[1];
sx q[1];
rz(-1.8073945) q[1];
sx q[1];
rz(0.7787849) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9433702) q[0];
sx q[0];
rz(-1.3226042) q[0];
sx q[0];
rz(-2.7944837) q[0];
rz(2.07431) q[2];
sx q[2];
rz(-0.52769606) q[2];
sx q[2];
rz(2.6380981) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(0.89214954) q[1];
sx q[1];
rz(-1.8037233) q[1];
sx q[1];
rz(1.1942785) q[1];
rz(-3.0770244) q[3];
sx q[3];
rz(-1.350292) q[3];
sx q[3];
rz(1.8622423) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-1.6918216) q[2];
sx q[2];
rz(-1.8372953) q[2];
sx q[2];
rz(2.4251535) q[2];
rz(1.6843494) q[3];
sx q[3];
rz(-1.4102035) q[3];
sx q[3];
rz(1.289207) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0338106) q[0];
sx q[0];
rz(-2.6066715) q[0];
sx q[0];
rz(0.15144908) q[0];
rz(2.9653213) q[1];
sx q[1];
rz(-1.9468032) q[1];
sx q[1];
rz(0.72296468) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3176214) q[0];
sx q[0];
rz(-0.184632) q[0];
sx q[0];
rz(-2.1872107) q[0];
rz(-pi) q[1];
x q[1];
rz(1.9480013) q[2];
sx q[2];
rz(-1.7282439) q[2];
sx q[2];
rz(2.8709656) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(2.0235698) q[1];
sx q[1];
rz(-1.6260864) q[1];
sx q[1];
rz(1.0667332) q[1];
rz(-pi) q[2];
rz(2.3530657) q[3];
sx q[3];
rz(-1.8619814) q[3];
sx q[3];
rz(2.9332719) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.67939776) q[2];
sx q[2];
rz(-1.858819) q[2];
sx q[2];
rz(-0.52552137) q[2];
rz(-2.8578791) q[3];
sx q[3];
rz(-1.107639) q[3];
sx q[3];
rz(0.39653683) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5491966) q[0];
sx q[0];
rz(-1.1239197) q[0];
sx q[0];
rz(-0.32604937) q[0];
rz(1.0992959) q[1];
sx q[1];
rz(-0.14930832) q[1];
sx q[1];
rz(2.2717448) q[1];
rz(-0.99156689) q[2];
sx q[2];
rz(-0.75800037) q[2];
sx q[2];
rz(-2.6108685) q[2];
rz(-0.59393926) q[3];
sx q[3];
rz(-0.38729061) q[3];
sx q[3];
rz(-0.95872986) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
