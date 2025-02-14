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
rz(0.99211168) q[0];
sx q[0];
rz(3.8881128) q[0];
sx q[0];
rz(9.9915656) q[0];
rz(1.2504638) q[1];
sx q[1];
rz(-1.992978) q[1];
sx q[1];
rz(2.221938) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.88094372) q[0];
sx q[0];
rz(-0.80505097) q[0];
sx q[0];
rz(2.9183003) q[0];
x q[1];
rz(-0.72534277) q[2];
sx q[2];
rz(-1.1769466) q[2];
sx q[2];
rz(0.79532901) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.39988841) q[1];
sx q[1];
rz(-2.2117858) q[1];
sx q[1];
rz(1.6287032) q[1];
rz(-1.7149431) q[3];
sx q[3];
rz(-1.0394154) q[3];
sx q[3];
rz(1.2856158) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(0.93996843) q[2];
sx q[2];
rz(-0.93279606) q[2];
sx q[2];
rz(2.0882108) q[2];
rz(-0.71422226) q[3];
sx q[3];
rz(-0.37319365) q[3];
sx q[3];
rz(1.8478954) q[3];
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
rz(-pi) q[0];
sx q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6222222) q[0];
sx q[0];
rz(-2.433233) q[0];
sx q[0];
rz(2.569662) q[0];
rz(1.2988623) q[1];
sx q[1];
rz(-2.3426901) q[1];
sx q[1];
rz(2.3518708) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8379254) q[0];
sx q[0];
rz(-1.2695489) q[0];
sx q[0];
rz(3.0409854) q[0];
rz(-pi) q[1];
rz(-1.4105807) q[2];
sx q[2];
rz(-0.92276697) q[2];
sx q[2];
rz(2.2367978) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-0.29282031) q[1];
sx q[1];
rz(-0.97859425) q[1];
sx q[1];
rz(1.2315537) q[1];
x q[2];
rz(0.67482194) q[3];
sx q[3];
rz(-1.9673507) q[3];
sx q[3];
rz(0.048585437) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(2.2907437) q[2];
sx q[2];
rz(-0.76883832) q[2];
sx q[2];
rz(-1.3624462) q[2];
rz(2.8507774) q[3];
sx q[3];
rz(-1.7751866) q[3];
sx q[3];
rz(3.0349558) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0027851) q[0];
sx q[0];
rz(-2.0464351) q[0];
sx q[0];
rz(2.441067) q[0];
rz(-2.6558212) q[1];
sx q[1];
rz(-2.3120717) q[1];
sx q[1];
rz(-2.9170759) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1248847) q[0];
sx q[0];
rz(-1.9210235) q[0];
sx q[0];
rz(-1.4136397) q[0];
x q[1];
rz(-1.7211368) q[2];
sx q[2];
rz(-2.5579961) q[2];
sx q[2];
rz(2.4395669) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-2.3697046) q[1];
sx q[1];
rz(-1.6460895) q[1];
sx q[1];
rz(2.1545707) q[1];
rz(-pi) q[2];
rz(1.443601) q[3];
sx q[3];
rz(-0.7326829) q[3];
sx q[3];
rz(-0.81058433) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(0.38480467) q[2];
sx q[2];
rz(-2.2497358) q[2];
sx q[2];
rz(3.0925114) q[2];
rz(-3.0745506) q[3];
sx q[3];
rz(-1.9837572) q[3];
sx q[3];
rz(-2.4750347) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9662358) q[0];
sx q[0];
rz(-0.55713621) q[0];
sx q[0];
rz(-0.019158451) q[0];
rz(-1.7823904) q[1];
sx q[1];
rz(-1.4298341) q[1];
sx q[1];
rz(3.137099) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2990103) q[0];
sx q[0];
rz(-1.0753462) q[0];
sx q[0];
rz(-2.5240077) q[0];
rz(-pi) q[1];
rz(1.121663) q[2];
sx q[2];
rz(-2.0259078) q[2];
sx q[2];
rz(3.1246076) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.7178065) q[1];
sx q[1];
rz(-0.39925925) q[1];
sx q[1];
rz(2.0171793) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.4779303) q[3];
sx q[3];
rz(-1.3680653) q[3];
sx q[3];
rz(2.7331238) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(2.6898474) q[2];
sx q[2];
rz(-1.0901901) q[2];
sx q[2];
rz(1.4534265) q[2];
rz(-1.2545741) q[3];
sx q[3];
rz(-1.6202319) q[3];
sx q[3];
rz(0.5622676) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[2];
x q[3];
rz(-pi/2) q[3];
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
rz(-0.8101863) q[0];
sx q[0];
rz(-1.9047381) q[0];
sx q[0];
rz(1.1619262) q[0];
rz(-0.76238531) q[1];
sx q[1];
rz(-2.1547909) q[1];
sx q[1];
rz(-0.42957482) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3439826) q[0];
sx q[0];
rz(-1.6071885) q[0];
sx q[0];
rz(-1.162536) q[0];
rz(-pi) q[1];
rz(1.8354906) q[2];
sx q[2];
rz(-2.3117723) q[2];
sx q[2];
rz(-0.13665567) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(0.030368806) q[1];
sx q[1];
rz(-2.3850523) q[1];
sx q[1];
rz(-1.0683505) q[1];
rz(-0.3394903) q[3];
sx q[3];
rz(-0.94565839) q[3];
sx q[3];
rz(-3.1200266) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(2.9261711) q[2];
sx q[2];
rz(-1.8502219) q[2];
sx q[2];
rz(-1.8298836) q[2];
rz(-0.46071509) q[3];
sx q[3];
rz(-2.1381133) q[3];
sx q[3];
rz(3.0084012) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
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
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8120414) q[0];
sx q[0];
rz(-2.9819745) q[0];
sx q[0];
rz(2.0352236) q[0];
rz(-0.1144935) q[1];
sx q[1];
rz(-1.0527) q[1];
sx q[1];
rz(-3.0879424) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0993529) q[0];
sx q[0];
rz(-1.5650378) q[0];
sx q[0];
rz(1.8264218) q[0];
rz(-pi) q[1];
x q[1];
rz(2.3340763) q[2];
sx q[2];
rz(-1.4984594) q[2];
sx q[2];
rz(-1.7010207) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.4380456) q[1];
sx q[1];
rz(-1.5168191) q[1];
sx q[1];
rz(-2.0825858) q[1];
x q[2];
rz(-0.76308454) q[3];
sx q[3];
rz(-0.79417919) q[3];
sx q[3];
rz(1.1556243) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-1.2758808) q[2];
sx q[2];
rz(-1.2578473) q[2];
sx q[2];
rz(2.6344521) q[2];
rz(1.3745314) q[3];
sx q[3];
rz(-1.5406698) q[3];
sx q[3];
rz(-1.9929632) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
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
rz(0.62992612) q[0];
sx q[0];
rz(-1.1277132) q[0];
sx q[0];
rz(0.041286904) q[0];
rz(2.4033026) q[1];
sx q[1];
rz(-1.9169044) q[1];
sx q[1];
rz(0.98947492) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.51136298) q[0];
sx q[0];
rz(-1.0697027) q[0];
sx q[0];
rz(2.5998235) q[0];
rz(-pi) q[1];
x q[1];
rz(0.47669741) q[2];
sx q[2];
rz(-0.39604353) q[2];
sx q[2];
rz(0.11644289) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-1.7478024) q[1];
sx q[1];
rz(-0.89266864) q[1];
sx q[1];
rz(0.80121471) q[1];
rz(-2.0416077) q[3];
sx q[3];
rz(-2.4306477) q[3];
sx q[3];
rz(2.8266852) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-1.5422633) q[2];
sx q[2];
rz(-1.0487391) q[2];
sx q[2];
rz(0.85912022) q[2];
rz(-3.1298992) q[3];
sx q[3];
rz(-1.6418567) q[3];
sx q[3];
rz(3.0288127) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.82666731) q[0];
sx q[0];
rz(-0.59190094) q[0];
sx q[0];
rz(2.9097606) q[0];
rz(-0.58206093) q[1];
sx q[1];
rz(-1.3287909) q[1];
sx q[1];
rz(1.2190762) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.40755475) q[0];
sx q[0];
rz(-2.1269607) q[0];
sx q[0];
rz(-1.4329877) q[0];
rz(-pi) q[1];
rz(-2.788576) q[2];
sx q[2];
rz(-1.5215877) q[2];
sx q[2];
rz(0.69248688) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-2.9207014) q[1];
sx q[1];
rz(-2.9013843) q[1];
sx q[1];
rz(2.2574053) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.3858504) q[3];
sx q[3];
rz(-1.1583987) q[3];
sx q[3];
rz(-0.88471593) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(1.1845188) q[2];
sx q[2];
rz(-1.3513214) q[2];
sx q[2];
rz(2.3547724) q[2];
rz(1.6400853) q[3];
sx q[3];
rz(-0.4314751) q[3];
sx q[3];
rz(1.4260346) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
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
rz(-2.1805434) q[0];
sx q[0];
rz(-1.3732055) q[0];
sx q[0];
rz(-2.2013262) q[0];
rz(2.6967948) q[1];
sx q[1];
rz(-2.2856789) q[1];
sx q[1];
rz(2.99517) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6165926) q[0];
sx q[0];
rz(-1.6967275) q[0];
sx q[0];
rz(0.30091826) q[0];
x q[1];
rz(1.8249885) q[2];
sx q[2];
rz(-1.3066402) q[2];
sx q[2];
rz(-3.11655) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-0.89288974) q[1];
sx q[1];
rz(-1.8161402) q[1];
sx q[1];
rz(-0.69458436) q[1];
x q[2];
rz(-1.4692471) q[3];
sx q[3];
rz(-0.78208215) q[3];
sx q[3];
rz(-2.8080432) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(3.0545097) q[2];
sx q[2];
rz(-1.4887709) q[2];
sx q[2];
rz(2.3040237) q[2];
rz(0.93291035) q[3];
sx q[3];
rz(-2.1541903) q[3];
sx q[3];
rz(0.78105175) q[3];
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
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8785716) q[0];
sx q[0];
rz(-1.913338) q[0];
sx q[0];
rz(-1.7735057) q[0];
rz(-0.076016501) q[1];
sx q[1];
rz(-0.58034211) q[1];
sx q[1];
rz(-1.1200294) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4809326) q[0];
sx q[0];
rz(-1.8802065) q[0];
sx q[0];
rz(-0.27746986) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.4330288) q[2];
sx q[2];
rz(-1.3363133) q[2];
sx q[2];
rz(-1.3402779) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(1.489913) q[1];
sx q[1];
rz(-1.8615684) q[1];
sx q[1];
rz(-3.0417697) q[1];
rz(-2.132776) q[3];
sx q[3];
rz(-0.19246092) q[3];
sx q[3];
rz(0.61103067) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(2.7398305) q[2];
sx q[2];
rz(-1.4515452) q[2];
sx q[2];
rz(0.24701992) q[2];
rz(-2.5714827) q[3];
sx q[3];
rz(-2.6542122) q[3];
sx q[3];
rz(2.0232239) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
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
rz(1.0982672) q[0];
sx q[0];
rz(-1.7541616) q[0];
sx q[0];
rz(2.7761205) q[0];
rz(-0.098516057) q[1];
sx q[1];
rz(-2.5010074) q[1];
sx q[1];
rz(0.88190257) q[1];
rz(-1.983704) q[2];
sx q[2];
rz(-1.558424) q[2];
sx q[2];
rz(1.3267714) q[2];
rz(2.5539342) q[3];
sx q[3];
rz(-0.44787221) q[3];
sx q[3];
rz(-0.032240562) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
