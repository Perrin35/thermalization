OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(0.87941909) q[0];
sx q[0];
rz(4.8882422) q[0];
sx q[0];
rz(12.56765) q[0];
rz(-1.6969504) q[1];
sx q[1];
rz(-2.0386219) q[1];
sx q[1];
rz(-2.3666518) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.66520663) q[0];
sx q[0];
rz(-1.5210266) q[0];
sx q[0];
rz(-0.14270356) q[0];
rz(-0.42595072) q[2];
sx q[2];
rz(-0.89086878) q[2];
sx q[2];
rz(1.8217063) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-0.98382681) q[1];
sx q[1];
rz(-2.9938934) q[1];
sx q[1];
rz(1.8450792) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.48304708) q[3];
sx q[3];
rz(-0.31937283) q[3];
sx q[3];
rz(-1.7628302) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(1.0455735) q[2];
sx q[2];
rz(-1.1117659) q[2];
sx q[2];
rz(-1.9457031) q[2];
rz(-1.1536417) q[3];
sx q[3];
rz(-2.3524645) q[3];
sx q[3];
rz(-1.7529863) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[2];
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
rz(-0.7213223) q[0];
sx q[0];
rz(-0.35636154) q[0];
sx q[0];
rz(1.8700245) q[0];
rz(-1.0999854) q[1];
sx q[1];
rz(-1.1106691) q[1];
sx q[1];
rz(1.3756479) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.401706) q[0];
sx q[0];
rz(-1.6580083) q[0];
sx q[0];
rz(2.7313822) q[0];
rz(-pi) q[1];
x q[1];
rz(0.045250821) q[2];
sx q[2];
rz(-1.3327193) q[2];
sx q[2];
rz(-2.0766052) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(1.1609636) q[1];
sx q[1];
rz(-1.3914319) q[1];
sx q[1];
rz(2.7705926) q[1];
rz(-pi) q[2];
x q[2];
rz(0.35956412) q[3];
sx q[3];
rz(-2.095795) q[3];
sx q[3];
rz(2.0294702) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-2.9282844) q[2];
sx q[2];
rz(-1.3457315) q[2];
sx q[2];
rz(2.6518872) q[2];
rz(-1.1335763) q[3];
sx q[3];
rz(-0.19530345) q[3];
sx q[3];
rz(1.0666696) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5415444) q[0];
sx q[0];
rz(-1.8627889) q[0];
sx q[0];
rz(0.24060732) q[0];
rz(-0.34257564) q[1];
sx q[1];
rz(-0.97476417) q[1];
sx q[1];
rz(-1.906357) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.19757195) q[0];
sx q[0];
rz(-2.4320142) q[0];
sx q[0];
rz(2.5463085) q[0];
rz(-pi) q[1];
rz(-1.3182993) q[2];
sx q[2];
rz(-1.2090948) q[2];
sx q[2];
rz(0.30644882) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(2.0551128) q[1];
sx q[1];
rz(-1.8243454) q[1];
sx q[1];
rz(0.26992814) q[1];
rz(-2.5127605) q[3];
sx q[3];
rz(-1.1662081) q[3];
sx q[3];
rz(-1.7189327) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-2.3774595) q[2];
sx q[2];
rz(-1.5185792) q[2];
sx q[2];
rz(-1.4366478) q[2];
rz(-1.7403587) q[3];
sx q[3];
rz(-1.2652206) q[3];
sx q[3];
rz(-0.60825545) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.075994611) q[0];
sx q[0];
rz(-1.5290715) q[0];
sx q[0];
rz(2.3377989) q[0];
rz(2.1919788) q[1];
sx q[1];
rz(-1.6751553) q[1];
sx q[1];
rz(3.0217357) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0474931) q[0];
sx q[0];
rz(-0.96512981) q[0];
sx q[0];
rz(-0.043461965) q[0];
rz(-1.4255964) q[2];
sx q[2];
rz(-1.2135266) q[2];
sx q[2];
rz(-1.6574588) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(2.4704628) q[1];
sx q[1];
rz(-1.5884591) q[1];
sx q[1];
rz(-1.3385685) q[1];
rz(-2.4265392) q[3];
sx q[3];
rz(-2.4063769) q[3];
sx q[3];
rz(-0.27763593) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(2.5084761) q[2];
sx q[2];
rz(-1.2458331) q[2];
sx q[2];
rz(-3.0692549) q[2];
rz(-2.7667601) q[3];
sx q[3];
rz(-2.4852677) q[3];
sx q[3];
rz(1.1981296) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4478093) q[0];
sx q[0];
rz(-2.5700975) q[0];
sx q[0];
rz(-2.6547292) q[0];
rz(2.411719) q[1];
sx q[1];
rz(-2.2327773) q[1];
sx q[1];
rz(-1.240085) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4362674) q[0];
sx q[0];
rz(-1.3316532) q[0];
sx q[0];
rz(1.5572085) q[0];
rz(-2.3104722) q[2];
sx q[2];
rz(-2.7925425) q[2];
sx q[2];
rz(2.8380053) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(2.005849) q[1];
sx q[1];
rz(-1.6372576) q[1];
sx q[1];
rz(0.60484109) q[1];
rz(2.2206743) q[3];
sx q[3];
rz(-0.71483597) q[3];
sx q[3];
rz(-0.20540796) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(2.038625) q[2];
sx q[2];
rz(-2.2323148) q[2];
sx q[2];
rz(-2.4385578) q[2];
rz(1.4098343) q[3];
sx q[3];
rz(-1.792428) q[3];
sx q[3];
rz(0.15771244) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.96419656) q[0];
sx q[0];
rz(-1.2921768) q[0];
sx q[0];
rz(-0.99037209) q[0];
rz(3.0888427) q[1];
sx q[1];
rz(-0.92664781) q[1];
sx q[1];
rz(1.6606768) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2213343) q[0];
sx q[0];
rz(-1.1898367) q[0];
sx q[0];
rz(-0.078698054) q[0];
x q[1];
rz(3.1408177) q[2];
sx q[2];
rz(-1.1256071) q[2];
sx q[2];
rz(-1.0926525) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(0.78696886) q[1];
sx q[1];
rz(-2.0957392) q[1];
sx q[1];
rz(2.900219) q[1];
rz(1.6310286) q[3];
sx q[3];
rz(-0.67959736) q[3];
sx q[3];
rz(0.70590245) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(0.0084373077) q[2];
sx q[2];
rz(-1.6494273) q[2];
sx q[2];
rz(-2.5793502) q[2];
rz(1.0605313) q[3];
sx q[3];
rz(-0.74354592) q[3];
sx q[3];
rz(0.26091584) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
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
rz(0.19787191) q[0];
sx q[0];
rz(-1.3523538) q[0];
sx q[0];
rz(1.4690171) q[0];
rz(-1.0143657) q[1];
sx q[1];
rz(-1.0275774) q[1];
sx q[1];
rz(-1.8168824) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6727407) q[0];
sx q[0];
rz(-1.6981089) q[0];
sx q[0];
rz(1.5479814) q[0];
rz(2.1340738) q[2];
sx q[2];
rz(-2.0118606) q[2];
sx q[2];
rz(1.2457459) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-1.3130256) q[1];
sx q[1];
rz(-1.1638068) q[1];
sx q[1];
rz(1.9810956) q[1];
x q[2];
rz(-2.7571194) q[3];
sx q[3];
rz(-0.95015804) q[3];
sx q[3];
rz(-1.4409161) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-1.5801195) q[2];
sx q[2];
rz(-1.7529528) q[2];
sx q[2];
rz(2.6980147) q[2];
rz(-2.1929072) q[3];
sx q[3];
rz(-1.9493999) q[3];
sx q[3];
rz(0.77478066) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
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
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.401944) q[0];
sx q[0];
rz(-2.1304603) q[0];
sx q[0];
rz(-2.6810714) q[0];
rz(0.1000239) q[1];
sx q[1];
rz(-0.99383751) q[1];
sx q[1];
rz(1.2896279) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.078482) q[0];
sx q[0];
rz(-1.0149628) q[0];
sx q[0];
rz(-2.1102935) q[0];
rz(2.3233534) q[2];
sx q[2];
rz(-2.1365676) q[2];
sx q[2];
rz(-0.73275369) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-2.928927) q[1];
sx q[1];
rz(-1.3839098) q[1];
sx q[1];
rz(3.0287663) q[1];
rz(-pi) q[2];
rz(1.4054221) q[3];
sx q[3];
rz(-0.91455063) q[3];
sx q[3];
rz(1.0429494) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-0.9926247) q[2];
sx q[2];
rz(-2.830539) q[2];
sx q[2];
rz(2.0588493) q[2];
rz(0.096171245) q[3];
sx q[3];
rz(-1.440719) q[3];
sx q[3];
rz(1.2307897) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7643395) q[0];
sx q[0];
rz(-2.9057673) q[0];
sx q[0];
rz(-0.33426958) q[0];
rz(-1.224068) q[1];
sx q[1];
rz(-1.5806438) q[1];
sx q[1];
rz(-0.28265488) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6982272) q[0];
sx q[0];
rz(-2.3400314) q[0];
sx q[0];
rz(2.482224) q[0];
rz(-pi) q[1];
x q[1];
rz(1.7868144) q[2];
sx q[2];
rz(-1.2470494) q[2];
sx q[2];
rz(-2.266989) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.64687318) q[1];
sx q[1];
rz(-0.44534007) q[1];
sx q[1];
rz(-0.41782197) q[1];
x q[2];
rz(1.4275527) q[3];
sx q[3];
rz(-2.4510265) q[3];
sx q[3];
rz(2.7018202) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(0.21215542) q[2];
sx q[2];
rz(-1.9781457) q[2];
sx q[2];
rz(0.7412509) q[2];
rz(-2.6397928) q[3];
sx q[3];
rz(-1.3391677) q[3];
sx q[3];
rz(-0.58399502) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.042645) q[0];
sx q[0];
rz(-2.2462923) q[0];
sx q[0];
rz(2.6328971) q[0];
rz(0.11518654) q[1];
sx q[1];
rz(-0.70786628) q[1];
sx q[1];
rz(-2.4597816) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6012321) q[0];
sx q[0];
rz(-1.0245748) q[0];
sx q[0];
rz(-0.19646125) q[0];
rz(-pi) q[1];
rz(2.290906) q[2];
sx q[2];
rz(-2.2518573) q[2];
sx q[2];
rz(-2.1249287) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(2.0360003) q[1];
sx q[1];
rz(-2.3899374) q[1];
sx q[1];
rz(0.73026258) q[1];
rz(-pi) q[2];
x q[2];
rz(0.52313157) q[3];
sx q[3];
rz(-2.0683859) q[3];
sx q[3];
rz(2.350972) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.8455785) q[2];
sx q[2];
rz(-1.6342376) q[2];
sx q[2];
rz(0.34109035) q[2];
rz(-2.0579445) q[3];
sx q[3];
rz(-0.79569474) q[3];
sx q[3];
rz(2.9705689) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[2];
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
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1682128) q[0];
sx q[0];
rz(-1.6642878) q[0];
sx q[0];
rz(-0.99933495) q[0];
rz(-2.534261) q[1];
sx q[1];
rz(-0.9098396) q[1];
sx q[1];
rz(0.20914016) q[1];
rz(0.63275679) q[2];
sx q[2];
rz(-0.7706332) q[2];
sx q[2];
rz(2.3494233) q[2];
rz(-1.7726462) q[3];
sx q[3];
rz(-1.1931843) q[3];
sx q[3];
rz(-1.9184792) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
