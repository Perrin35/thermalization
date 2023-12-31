OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(1.4484654) q[0];
sx q[0];
rz(-2.6187596) q[0];
sx q[0];
rz(-2.5180106) q[0];
rz(-2.8514255) q[1];
sx q[1];
rz(-0.71915141) q[1];
sx q[1];
rz(-2.6410988) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.96007632) q[0];
sx q[0];
rz(-1.5746563) q[0];
sx q[0];
rz(0.5998248) q[0];
x q[1];
rz(-1.4209461) q[2];
sx q[2];
rz(-1.4784001) q[2];
sx q[2];
rz(2.7211651) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(3.0296368) q[1];
sx q[1];
rz(-1.2458548) q[1];
sx q[1];
rz(-0.8128266) q[1];
rz(1.2455363) q[3];
sx q[3];
rz(-0.45913011) q[3];
sx q[3];
rz(-1.1577275) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-1.3502675) q[2];
sx q[2];
rz(-1.2056377) q[2];
sx q[2];
rz(1.9187437) q[2];
rz(1.4482927) q[3];
sx q[3];
rz(-2.1494614) q[3];
sx q[3];
rz(-2.1712415) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7704849) q[0];
sx q[0];
rz(-1.4828232) q[0];
sx q[0];
rz(2.0626542) q[0];
rz(-1.7547912) q[1];
sx q[1];
rz(-2.3290122) q[1];
sx q[1];
rz(2.4761377) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5011713) q[0];
sx q[0];
rz(-1.2873532) q[0];
sx q[0];
rz(2.6675176) q[0];
rz(0.47768728) q[2];
sx q[2];
rz(-2.7825232) q[2];
sx q[2];
rz(-1.2815086) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(2.7484819) q[1];
sx q[1];
rz(-2.0512274) q[1];
sx q[1];
rz(-1.3039116) q[1];
rz(-pi) q[2];
rz(-2.2802248) q[3];
sx q[3];
rz(-1.5449761) q[3];
sx q[3];
rz(0.65755075) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(0.60454303) q[2];
sx q[2];
rz(-1.4031354) q[2];
sx q[2];
rz(-0.075142168) q[2];
rz(-1.6710619) q[3];
sx q[3];
rz(-2.0139147) q[3];
sx q[3];
rz(2.4501734) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.55494088) q[0];
sx q[0];
rz(-1.8692769) q[0];
sx q[0];
rz(-0.75876045) q[0];
rz(1.2930019) q[1];
sx q[1];
rz(-1.2932152) q[1];
sx q[1];
rz(-1.4000777) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5132719) q[0];
sx q[0];
rz(-1.3043881) q[0];
sx q[0];
rz(0.4011641) q[0];
rz(-pi) q[1];
x q[1];
rz(2.5163469) q[2];
sx q[2];
rz(-2.9125104) q[2];
sx q[2];
rz(-2.2428227) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(1.2565508) q[1];
sx q[1];
rz(-0.98837822) q[1];
sx q[1];
rz(0.38862733) q[1];
rz(-pi) q[2];
x q[2];
rz(0.16617822) q[3];
sx q[3];
rz(-1.2697392) q[3];
sx q[3];
rz(0.64134669) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(2.489958) q[2];
sx q[2];
rz(-0.48214665) q[2];
sx q[2];
rz(0.65762323) q[2];
rz(1.970132) q[3];
sx q[3];
rz(-1.6572584) q[3];
sx q[3];
rz(1.7224147) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3477429) q[0];
sx q[0];
rz(-0.94809735) q[0];
sx q[0];
rz(1.4452274) q[0];
rz(-1.4472648) q[1];
sx q[1];
rz(-1.4935962) q[1];
sx q[1];
rz(0.34805527) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.30663438) q[0];
sx q[0];
rz(-1.8638896) q[0];
sx q[0];
rz(-2.4819863) q[0];
rz(2.9910827) q[2];
sx q[2];
rz(-1.2005271) q[2];
sx q[2];
rz(-2.9589257) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(3.0192249) q[1];
sx q[1];
rz(-2.5248563) q[1];
sx q[1];
rz(1.7255746) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.4481941) q[3];
sx q[3];
rz(-0.20892538) q[3];
sx q[3];
rz(-2.8914176) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(0.41670123) q[2];
sx q[2];
rz(-1.3092224) q[2];
sx q[2];
rz(0.42281881) q[2];
rz(0.73741284) q[3];
sx q[3];
rz(-0.80602065) q[3];
sx q[3];
rz(0.084658682) q[3];
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
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2622862) q[0];
sx q[0];
rz(-1.3129741) q[0];
sx q[0];
rz(-0.53043956) q[0];
rz(-2.2166705) q[1];
sx q[1];
rz(-1.5735807) q[1];
sx q[1];
rz(1.8431429) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0153774) q[0];
sx q[0];
rz(-1.4091638) q[0];
sx q[0];
rz(0.14590185) q[0];
rz(1.2735882) q[2];
sx q[2];
rz(-0.98515918) q[2];
sx q[2];
rz(2.1538018) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-3.0188705) q[1];
sx q[1];
rz(-1.4713305) q[1];
sx q[1];
rz(2.781267) q[1];
rz(0.88921806) q[3];
sx q[3];
rz(-2.4272356) q[3];
sx q[3];
rz(2.0444972) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(1.2003145) q[2];
sx q[2];
rz(-1.986074) q[2];
sx q[2];
rz(0.47362622) q[2];
rz(-3.04223) q[3];
sx q[3];
rz(-1.8615581) q[3];
sx q[3];
rz(0.84053269) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6376003) q[0];
sx q[0];
rz(-1.3562599) q[0];
sx q[0];
rz(-0.9978869) q[0];
rz(-2.2672794) q[1];
sx q[1];
rz(-2.120178) q[1];
sx q[1];
rz(-2.6748437) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.1186819) q[0];
sx q[0];
rz(-0.70436275) q[0];
sx q[0];
rz(-0.33606152) q[0];
rz(-pi) q[1];
rz(-2.7141063) q[2];
sx q[2];
rz(-1.7726521) q[2];
sx q[2];
rz(-2.0919378) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-3.0007243) q[1];
sx q[1];
rz(-1.6182401) q[1];
sx q[1];
rz(2.9904757) q[1];
rz(-pi) q[2];
x q[2];
rz(1.407133) q[3];
sx q[3];
rz(-1.4758849) q[3];
sx q[3];
rz(-1.519219) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-2.2816887) q[2];
sx q[2];
rz(-2.6624661) q[2];
sx q[2];
rz(1.5768645) q[2];
rz(2.5148897) q[3];
sx q[3];
rz(-1.3650711) q[3];
sx q[3];
rz(-0.68157649) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
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
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.80207434) q[0];
sx q[0];
rz(-3.0537362) q[0];
sx q[0];
rz(3.0457892) q[0];
rz(-pi) q[1];
rz(2.1599342) q[2];
sx q[2];
rz(-0.33827153) q[2];
sx q[2];
rz(2.8209518) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.35510264) q[1];
sx q[1];
rz(-0.10663248) q[1];
sx q[1];
rz(2.998888) q[1];
rz(0.94957955) q[3];
sx q[3];
rz(-2.5887244) q[3];
sx q[3];
rz(2.051193) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-0.55591136) q[2];
sx q[2];
rz(-2.6101117) q[2];
sx q[2];
rz(-1.4038203) q[2];
rz(0.79706556) q[3];
sx q[3];
rz(-2.6795487) q[3];
sx q[3];
rz(1.9246624) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7109011) q[0];
sx q[0];
rz(-0.71390188) q[0];
sx q[0];
rz(-0.28924334) q[0];
rz(-2.5166683) q[1];
sx q[1];
rz(-1.9754675) q[1];
sx q[1];
rz(1.3141059) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4411366) q[0];
sx q[0];
rz(-1.4608129) q[0];
sx q[0];
rz(0.76018795) q[0];
rz(-1.9004702) q[2];
sx q[2];
rz(-2.2909819) q[2];
sx q[2];
rz(-0.61851293) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.9574979) q[1];
sx q[1];
rz(-3.0093319) q[1];
sx q[1];
rz(1.1485419) q[1];
x q[2];
rz(-0.81359158) q[3];
sx q[3];
rz(-1.3435257) q[3];
sx q[3];
rz(0.0031301216) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(0.73359314) q[2];
sx q[2];
rz(-1.5780129) q[2];
sx q[2];
rz(-1.4286263) q[2];
rz(0.96380487) q[3];
sx q[3];
rz(-2.0644085) q[3];
sx q[3];
rz(2.111964) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.52255094) q[0];
sx q[0];
rz(-1.4551117) q[0];
sx q[0];
rz(-1.2458941) q[0];
rz(3.030581) q[1];
sx q[1];
rz(-1.9440034) q[1];
sx q[1];
rz(-2.5949809) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.36771691) q[0];
sx q[0];
rz(-1.9426632) q[0];
sx q[0];
rz(2.0433776) q[0];
x q[1];
rz(-1.941628) q[2];
sx q[2];
rz(-1.9413345) q[2];
sx q[2];
rz(1.759699) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-2.6555772) q[1];
sx q[1];
rz(-1.7421107) q[1];
sx q[1];
rz(-1.1737215) q[1];
rz(-pi) q[2];
rz(0.48874493) q[3];
sx q[3];
rz(-2.1663323) q[3];
sx q[3];
rz(0.83317703) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(2.1116011) q[2];
sx q[2];
rz(-2.2932055) q[2];
sx q[2];
rz(1.6938422) q[2];
rz(-2.0041806) q[3];
sx q[3];
rz(-1.3237938) q[3];
sx q[3];
rz(-0.64731961) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2820213) q[0];
sx q[0];
rz(-2.8347926) q[0];
sx q[0];
rz(2.4243673) q[0];
rz(-1.2099129) q[1];
sx q[1];
rz(-0.33214339) q[1];
sx q[1];
rz(2.4338914) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.087698547) q[0];
sx q[0];
rz(-2.2367034) q[0];
sx q[0];
rz(0.44184394) q[0];
rz(-pi) q[1];
x q[1];
rz(1.3466481) q[2];
sx q[2];
rz(-1.7094311) q[2];
sx q[2];
rz(-0.56425205) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-1.3860491) q[1];
sx q[1];
rz(-1.886133) q[1];
sx q[1];
rz(0.83868933) q[1];
x q[2];
rz(0.95146146) q[3];
sx q[3];
rz(-1.8639495) q[3];
sx q[3];
rz(-2.4952863) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(1.5132961) q[2];
sx q[2];
rz(-0.6568903) q[2];
sx q[2];
rz(-0.90325242) q[2];
rz(1.5385657) q[3];
sx q[3];
rz(-2.2730946) q[3];
sx q[3];
rz(-2.2911151) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.83508867) q[0];
sx q[0];
rz(-0.36515129) q[0];
sx q[0];
rz(-0.93602244) q[0];
rz(-0.8159591) q[1];
sx q[1];
rz(-0.42146704) q[1];
sx q[1];
rz(-2.0889919) q[1];
rz(1.1076526) q[2];
sx q[2];
rz(-1.9911498) q[2];
sx q[2];
rz(-0.1291612) q[2];
rz(0.64745263) q[3];
sx q[3];
rz(-3.0734607) q[3];
sx q[3];
rz(-1.4510696) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
