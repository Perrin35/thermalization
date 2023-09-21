OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(2.6089132) q[0];
sx q[0];
rz(-0.37663868) q[0];
sx q[0];
rz(0.11178804) q[0];
rz(1.6821661) q[1];
sx q[1];
rz(4.7987727) q[1];
sx q[1];
rz(6.12943) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.35225866) q[0];
sx q[0];
rz(-0.52868045) q[0];
sx q[0];
rz(2.5369011) q[0];
rz(-1.8122187) q[2];
sx q[2];
rz(-1.4960939) q[2];
sx q[2];
rz(-2.9818997) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-2.0220713) q[1];
sx q[1];
rz(-1.6011366) q[1];
sx q[1];
rz(-1.2519217) q[1];
x q[2];
rz(-2.2114803) q[3];
sx q[3];
rz(-2.6112587) q[3];
sx q[3];
rz(1.9213284) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-1.443632) q[2];
sx q[2];
rz(-1.4322832) q[2];
sx q[2];
rz(-1.704818) q[2];
rz(-0.73389655) q[3];
sx q[3];
rz(-1.5926444) q[3];
sx q[3];
rz(0.51600391) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2519418) q[0];
sx q[0];
rz(-1.2263068) q[0];
sx q[0];
rz(2.2170128) q[0];
rz(-0.997116) q[1];
sx q[1];
rz(-0.50874248) q[1];
sx q[1];
rz(1.8181713) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.91081496) q[0];
sx q[0];
rz(-1.9998904) q[0];
sx q[0];
rz(2.6399415) q[0];
rz(-pi) q[1];
x q[1];
rz(0.1214059) q[2];
sx q[2];
rz(-2.7896023) q[2];
sx q[2];
rz(-0.44167232) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-1.1393226) q[1];
sx q[1];
rz(-1.3011258) q[1];
sx q[1];
rz(-2.3766999) q[1];
rz(2.7534361) q[3];
sx q[3];
rz(-2.1106488) q[3];
sx q[3];
rz(-2.848958) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.20415846) q[2];
sx q[2];
rz(-1.5519451) q[2];
sx q[2];
rz(-2.3834174) q[2];
rz(-0.6289064) q[3];
sx q[3];
rz(-2.7401676) q[3];
sx q[3];
rz(-1.1531856) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4085098) q[0];
sx q[0];
rz(-1.874431) q[0];
sx q[0];
rz(-0.68840233) q[0];
rz(-0.06772659) q[1];
sx q[1];
rz(-1.3893145) q[1];
sx q[1];
rz(-2.6115131) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.63307525) q[0];
sx q[0];
rz(-1.2808262) q[0];
sx q[0];
rz(-3.0301208) q[0];
rz(1.0447787) q[2];
sx q[2];
rz(-1.9472497) q[2];
sx q[2];
rz(1.2029755) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-0.59280076) q[1];
sx q[1];
rz(-0.93344102) q[1];
sx q[1];
rz(2.7194276) q[1];
x q[2];
rz(-1.8102874) q[3];
sx q[3];
rz(-1.665984) q[3];
sx q[3];
rz(-2.5483607) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-2.3337341) q[2];
sx q[2];
rz(-3.1299751) q[2];
sx q[2];
rz(0.90144908) q[2];
rz(0.83550134) q[3];
sx q[3];
rz(-1.6136026) q[3];
sx q[3];
rz(-1.8301331) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
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
rz(2.9505342) q[0];
sx q[0];
rz(-1.5427417) q[0];
sx q[0];
rz(0.78432551) q[0];
rz(0.061231881) q[1];
sx q[1];
rz(-2.4274554) q[1];
sx q[1];
rz(3.004946) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.89350677) q[0];
sx q[0];
rz(-2.6419905) q[0];
sx q[0];
rz(1.8462371) q[0];
rz(-pi) q[1];
x q[1];
rz(2.6776671) q[2];
sx q[2];
rz(-1.8564965) q[2];
sx q[2];
rz(0.49109101) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(0.54388753) q[1];
sx q[1];
rz(-1.6391616) q[1];
sx q[1];
rz(-0.023936546) q[1];
rz(-pi) q[2];
rz(-2.1468048) q[3];
sx q[3];
rz(-1.4928865) q[3];
sx q[3];
rz(-0.60914492) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(2.1293929) q[2];
sx q[2];
rz(-2.2012074) q[2];
sx q[2];
rz(-2.5811035) q[2];
rz(-3.1292606) q[3];
sx q[3];
rz(-0.90356946) q[3];
sx q[3];
rz(-1.0906609) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
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
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7085768) q[0];
sx q[0];
rz(-0.60537678) q[0];
sx q[0];
rz(-0.82114712) q[0];
rz(-2.2654146) q[1];
sx q[1];
rz(-2.2416302) q[1];
sx q[1];
rz(-1.7339773) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6038937) q[0];
sx q[0];
rz(-1.1802117) q[0];
sx q[0];
rz(1.9498528) q[0];
rz(2.1528835) q[2];
sx q[2];
rz(-1.8277797) q[2];
sx q[2];
rz(0.35494057) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-2.447532) q[1];
sx q[1];
rz(-0.58528712) q[1];
sx q[1];
rz(1.8156169) q[1];
rz(-pi) q[2];
rz(-2.1122583) q[3];
sx q[3];
rz(-2.7280305) q[3];
sx q[3];
rz(-3.001861) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-1.6901107) q[2];
sx q[2];
rz(-1.9268945) q[2];
sx q[2];
rz(3.0991128) q[2];
rz(2.5111607) q[3];
sx q[3];
rz(-0.63215956) q[3];
sx q[3];
rz(0.49155864) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7522488) q[0];
sx q[0];
rz(-1.9135973) q[0];
sx q[0];
rz(2.65843) q[0];
rz(-1.0522316) q[1];
sx q[1];
rz(-1.9960884) q[1];
sx q[1];
rz(0.56484708) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7171377) q[0];
sx q[0];
rz(-2.0949445) q[0];
sx q[0];
rz(-0.15630396) q[0];
x q[1];
rz(1.1499314) q[2];
sx q[2];
rz(-1.5016342) q[2];
sx q[2];
rz(1.7771306) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(0.55826742) q[1];
sx q[1];
rz(-0.43154432) q[1];
sx q[1];
rz(-1.4285018) q[1];
x q[2];
rz(-0.23128831) q[3];
sx q[3];
rz(-1.7219208) q[3];
sx q[3];
rz(-0.95201991) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-0.57006449) q[2];
sx q[2];
rz(-1.0674942) q[2];
sx q[2];
rz(1.8035536) q[2];
rz(-1.3048874) q[3];
sx q[3];
rz(-1.0675425) q[3];
sx q[3];
rz(-0.6750955) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
x q[2];
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
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.17334443) q[0];
sx q[0];
rz(-1.4302379) q[0];
sx q[0];
rz(-2.5937953) q[0];
rz(-2.3563747) q[1];
sx q[1];
rz(-1.806587) q[1];
sx q[1];
rz(0.26842591) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1584394) q[0];
sx q[0];
rz(-1.8162677) q[0];
sx q[0];
rz(1.655274) q[0];
x q[1];
rz(0.78598778) q[2];
sx q[2];
rz(-0.95324264) q[2];
sx q[2];
rz(2.3812889) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-0.35560289) q[1];
sx q[1];
rz(-0.59521788) q[1];
sx q[1];
rz(1.7463513) q[1];
rz(1.0309585) q[3];
sx q[3];
rz(-1.3098048) q[3];
sx q[3];
rz(0.48418448) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(0.61775529) q[2];
sx q[2];
rz(-2.338151) q[2];
sx q[2];
rz(0.8141554) q[2];
rz(2.7653149) q[3];
sx q[3];
rz(-1.1637996) q[3];
sx q[3];
rz(3.0686839) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
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
rz(-2.5557264) q[0];
sx q[0];
rz(-1.4050452) q[0];
sx q[0];
rz(-2.1222173) q[0];
rz(-2.2881919) q[1];
sx q[1];
rz(-1.1420206) q[1];
sx q[1];
rz(-0.44874915) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.42496029) q[0];
sx q[0];
rz(-0.53629959) q[0];
sx q[0];
rz(2.0400356) q[0];
rz(0.25787392) q[2];
sx q[2];
rz(-1.8631301) q[2];
sx q[2];
rz(2.984798) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-1.6905578) q[1];
sx q[1];
rz(-2.096855) q[1];
sx q[1];
rz(-1.2760389) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.429871) q[3];
sx q[3];
rz(-2.2300643) q[3];
sx q[3];
rz(-0.014051138) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-2.2723508) q[2];
sx q[2];
rz(-1.7607471) q[2];
sx q[2];
rz(-2.8273919) q[2];
rz(-0.82434404) q[3];
sx q[3];
rz(-2.6894675) q[3];
sx q[3];
rz(-0.79469386) q[3];
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
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.070351275) q[0];
sx q[0];
rz(-3.0817139) q[0];
sx q[0];
rz(-1.2605793) q[0];
rz(2.4977327) q[1];
sx q[1];
rz(-1.9088129) q[1];
sx q[1];
rz(0.018928122) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8277055) q[0];
sx q[0];
rz(-2.8773617) q[0];
sx q[0];
rz(-0.068473579) q[0];
rz(-pi) q[1];
rz(1.6690427) q[2];
sx q[2];
rz(-1.4155404) q[2];
sx q[2];
rz(-1.9948024) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(0.70871204) q[1];
sx q[1];
rz(-1.6212665) q[1];
sx q[1];
rz(-2.9570079) q[1];
rz(1.8065679) q[3];
sx q[3];
rz(-1.8528432) q[3];
sx q[3];
rz(0.75197938) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-0.54946047) q[2];
sx q[2];
rz(-0.38828725) q[2];
sx q[2];
rz(2.4712759) q[2];
rz(-2.629225) q[3];
sx q[3];
rz(-1.7497601) q[3];
sx q[3];
rz(1.4204773) q[3];
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
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5877514) q[0];
sx q[0];
rz(-1.2441664) q[0];
sx q[0];
rz(-2.642139) q[0];
rz(-1.5746501) q[1];
sx q[1];
rz(-2.8630239) q[1];
sx q[1];
rz(2.0589028) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.85031803) q[0];
sx q[0];
rz(-0.58915888) q[0];
sx q[0];
rz(0.088081443) q[0];
x q[1];
rz(-0.81929368) q[2];
sx q[2];
rz(-1.5466006) q[2];
sx q[2];
rz(1.0011315) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-0.078552695) q[1];
sx q[1];
rz(-1.5895956) q[1];
sx q[1];
rz(2.5979554) q[1];
x q[2];
rz(-0.70268481) q[3];
sx q[3];
rz(-0.48326884) q[3];
sx q[3];
rz(2.0600704) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-1.3988951) q[2];
sx q[2];
rz(-1.164914) q[2];
sx q[2];
rz(2.6297074) q[2];
rz(0.39294696) q[3];
sx q[3];
rz(-1.4018551) q[3];
sx q[3];
rz(-2.0846562) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.95505161) q[0];
sx q[0];
rz(-1.3165836) q[0];
sx q[0];
rz(-2.5008428) q[0];
rz(2.4004249) q[1];
sx q[1];
rz(-2.3186431) q[1];
sx q[1];
rz(2.9021312) q[1];
rz(2.7411186) q[2];
sx q[2];
rz(-2.1486712) q[2];
sx q[2];
rz(1.490544) q[2];
rz(2.3951204) q[3];
sx q[3];
rz(-1.6975879) q[3];
sx q[3];
rz(-3.0293037) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];