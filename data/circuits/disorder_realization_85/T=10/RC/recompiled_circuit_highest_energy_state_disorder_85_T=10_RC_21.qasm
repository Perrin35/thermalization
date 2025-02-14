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
rz(-2.149481) q[0];
sx q[0];
rz(-0.7465201) q[0];
sx q[0];
rz(-0.56678766) q[0];
rz(-1.8911288) q[1];
sx q[1];
rz(-1.1486147) q[1];
sx q[1];
rz(-2.221938) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5772668) q[0];
sx q[0];
rz(-2.3503605) q[0];
sx q[0];
rz(-1.7971695) q[0];
rz(2.5819725) q[2];
sx q[2];
rz(-2.3336448) q[2];
sx q[2];
rz(-1.9576278) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.2055605) q[1];
sx q[1];
rz(-1.6171997) q[1];
sx q[1];
rz(0.64179365) q[1];
x q[2];
rz(0.53594671) q[3];
sx q[3];
rz(-1.4466373) q[3];
sx q[3];
rz(-0.21175993) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(0.93996843) q[2];
sx q[2];
rz(-0.93279606) q[2];
sx q[2];
rz(-1.0533818) q[2];
rz(-0.71422226) q[3];
sx q[3];
rz(-2.768399) q[3];
sx q[3];
rz(-1.8478954) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6222222) q[0];
sx q[0];
rz(-2.433233) q[0];
sx q[0];
rz(0.57193065) q[0];
rz(-1.2988623) q[1];
sx q[1];
rz(-0.79890257) q[1];
sx q[1];
rz(-0.78972185) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.1173218) q[0];
sx q[0];
rz(-0.3171176) q[0];
sx q[0];
rz(1.2581403) q[0];
x q[1];
rz(2.9339173) q[2];
sx q[2];
rz(-0.66476124) q[2];
sx q[2];
rz(-1.1663811) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-0.27093664) q[1];
sx q[1];
rz(-2.4693477) q[1];
sx q[1];
rz(-2.6822374) q[1];
rz(-pi) q[2];
x q[2];
rz(2.5511207) q[3];
sx q[3];
rz(-2.3749224) q[3];
sx q[3];
rz(-1.1696512) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(0.85084891) q[2];
sx q[2];
rz(-0.76883832) q[2];
sx q[2];
rz(1.3624462) q[2];
rz(0.29081523) q[3];
sx q[3];
rz(-1.7751866) q[3];
sx q[3];
rz(0.10663685) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
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
rz(-3.0027851) q[0];
sx q[0];
rz(-2.0464351) q[0];
sx q[0];
rz(-2.441067) q[0];
rz(-0.48577148) q[1];
sx q[1];
rz(-2.3120717) q[1];
sx q[1];
rz(-0.22451678) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.49977127) q[0];
sx q[0];
rz(-1.4232521) q[0];
sx q[0];
rz(0.35424252) q[0];
rz(1.7211368) q[2];
sx q[2];
rz(-0.58359658) q[2];
sx q[2];
rz(-0.70202578) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-2.3697046) q[1];
sx q[1];
rz(-1.4955031) q[1];
sx q[1];
rz(-2.1545707) q[1];
x q[2];
rz(-3.0279474) q[3];
sx q[3];
rz(-0.84535852) q[3];
sx q[3];
rz(2.501373) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(2.756788) q[2];
sx q[2];
rz(-2.2497358) q[2];
sx q[2];
rz(0.049081238) q[2];
rz(-3.0745506) q[3];
sx q[3];
rz(-1.9837572) q[3];
sx q[3];
rz(0.66655794) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1753569) q[0];
sx q[0];
rz(-0.55713621) q[0];
sx q[0];
rz(-3.1224342) q[0];
rz(1.7823904) q[1];
sx q[1];
rz(-1.4298341) q[1];
sx q[1];
rz(0.0044936831) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.27994) q[0];
sx q[0];
rz(-2.3706382) q[0];
sx q[0];
rz(-0.75087913) q[0];
rz(-pi) q[1];
rz(2.0199297) q[2];
sx q[2];
rz(-2.0259078) q[2];
sx q[2];
rz(-3.1246076) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(1.5792721) q[1];
sx q[1];
rz(-1.4021789) q[1];
sx q[1];
rz(1.9344485) q[1];
rz(-pi) q[2];
x q[2];
rz(0.2035844) q[3];
sx q[3];
rz(-1.4798375) q[3];
sx q[3];
rz(-1.1435777) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.4517453) q[2];
sx q[2];
rz(-1.0901901) q[2];
sx q[2];
rz(-1.6881662) q[2];
rz(-1.8870185) q[3];
sx q[3];
rz(-1.6202319) q[3];
sx q[3];
rz(2.5793251) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
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
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.8101863) q[0];
sx q[0];
rz(-1.2368546) q[0];
sx q[0];
rz(-1.9796665) q[0];
rz(0.76238531) q[1];
sx q[1];
rz(-0.98680174) q[1];
sx q[1];
rz(2.7120178) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3439826) q[0];
sx q[0];
rz(-1.6071885) q[0];
sx q[0];
rz(-1.162536) q[0];
x q[1];
rz(-7/(8*pi)) q[2];
sx q[2];
rz(-0.77829153) q[2];
sx q[2];
rz(-0.24519224) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(1.9813198) q[1];
sx q[1];
rz(-1.2339051) q[1];
sx q[1];
rz(0.87967061) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.1384835) q[3];
sx q[3];
rz(-2.4412182) q[3];
sx q[3];
rz(-0.56454235) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-2.9261711) q[2];
sx q[2];
rz(-1.8502219) q[2];
sx q[2];
rz(1.8298836) q[2];
rz(-2.6808776) q[3];
sx q[3];
rz(-1.0034794) q[3];
sx q[3];
rz(-0.13319143) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.32955125) q[0];
sx q[0];
rz(-0.1596182) q[0];
sx q[0];
rz(1.106369) q[0];
rz(0.1144935) q[1];
sx q[1];
rz(-2.0888927) q[1];
sx q[1];
rz(0.053650275) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.44941266) q[0];
sx q[0];
rz(-2.8859038) q[0];
sx q[0];
rz(-1.5480255) q[0];
rz(-pi) q[1];
x q[1];
rz(-3.0416415) q[2];
sx q[2];
rz(-2.3315773) q[2];
sx q[2];
rz(-0.19925995) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-0.96300321) q[1];
sx q[1];
rz(-0.51437639) q[1];
sx q[1];
rz(1.4609171) q[1];
rz(-0.63416449) q[3];
sx q[3];
rz(-2.0863219) q[3];
sx q[3];
rz(-0.17532119) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(1.8657118) q[2];
sx q[2];
rz(-1.2578473) q[2];
sx q[2];
rz(-0.50714058) q[2];
rz(1.7670613) q[3];
sx q[3];
rz(-1.5406698) q[3];
sx q[3];
rz(1.9929632) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.62992612) q[0];
sx q[0];
rz(-1.1277132) q[0];
sx q[0];
rz(-3.1003057) q[0];
rz(-0.73829007) q[1];
sx q[1];
rz(-1.2246882) q[1];
sx q[1];
rz(-0.98947492) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.51136298) q[0];
sx q[0];
rz(-1.0697027) q[0];
sx q[0];
rz(-2.5998235) q[0];
rz(-pi) q[1];
rz(2.7858753) q[2];
sx q[2];
rz(-1.748744) q[2];
sx q[2];
rz(-1.8989814) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-1.7478024) q[1];
sx q[1];
rz(-2.248924) q[1];
sx q[1];
rz(2.3403779) q[1];
rz(-pi) q[2];
rz(-2.0416077) q[3];
sx q[3];
rz(-0.71094497) q[3];
sx q[3];
rz(0.31490745) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-1.5993293) q[2];
sx q[2];
rz(-2.0928536) q[2];
sx q[2];
rz(-2.2824724) q[2];
rz(3.1298992) q[3];
sx q[3];
rz(-1.499736) q[3];
sx q[3];
rz(3.0288127) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.82666731) q[0];
sx q[0];
rz(-2.5496917) q[0];
sx q[0];
rz(2.9097606) q[0];
rz(2.5595317) q[1];
sx q[1];
rz(-1.3287909) q[1];
sx q[1];
rz(1.2190762) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9052637) q[0];
sx q[0];
rz(-1.6877315) q[0];
sx q[0];
rz(-0.56044436) q[0];
rz(-pi) q[1];
rz(-0.14149551) q[2];
sx q[2];
rz(-0.35628755) q[2];
sx q[2];
rz(-2.3959999) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(2.2197709) q[1];
sx q[1];
rz(-1.3857462) q[1];
sx q[1];
rz(2.9875523) q[1];
rz(-pi) q[2];
x q[2];
rz(0.7557423) q[3];
sx q[3];
rz(-1.983194) q[3];
sx q[3];
rz(0.88471593) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(1.1845188) q[2];
sx q[2];
rz(-1.3513214) q[2];
sx q[2];
rz(-0.78682023) q[2];
rz(-1.6400853) q[3];
sx q[3];
rz(-2.7101176) q[3];
sx q[3];
rz(-1.7155581) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1805434) q[0];
sx q[0];
rz(-1.7683872) q[0];
sx q[0];
rz(2.2013262) q[0];
rz(2.6967948) q[1];
sx q[1];
rz(-0.85591379) q[1];
sx q[1];
rz(-2.99517) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.43035591) q[0];
sx q[0];
rz(-0.32546639) q[0];
sx q[0];
rz(-0.40367608) q[0];
x q[1];
rz(-2.3926061) q[2];
sx q[2];
rz(-2.7770713) q[2];
sx q[2];
rz(-0.75789343) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(2.1799911) q[1];
sx q[1];
rz(-2.4117984) q[1];
sx q[1];
rz(0.37288937) q[1];
x q[2];
rz(-0.10036631) q[3];
sx q[3];
rz(-2.3477738) q[3];
sx q[3];
rz(-0.19089261) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-3.0545097) q[2];
sx q[2];
rz(-1.4887709) q[2];
sx q[2];
rz(-2.3040237) q[2];
rz(-0.93291035) q[3];
sx q[3];
rz(-2.1541903) q[3];
sx q[3];
rz(2.3605409) q[3];
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
sx q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.26302108) q[0];
sx q[0];
rz(-1.913338) q[0];
sx q[0];
rz(1.7735057) q[0];
rz(3.0655762) q[1];
sx q[1];
rz(-0.58034211) q[1];
sx q[1];
rz(-1.1200294) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4809326) q[0];
sx q[0];
rz(-1.8802065) q[0];
sx q[0];
rz(0.27746986) q[0];
x q[1];
rz(-0.35181184) q[2];
sx q[2];
rz(-0.73992554) q[2];
sx q[2];
rz(-0.49525317) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-1.489913) q[1];
sx q[1];
rz(-1.2800242) q[1];
sx q[1];
rz(-3.0417697) q[1];
rz(-pi) q[2];
rz(-3.038123) q[3];
sx q[3];
rz(-1.7333687) q[3];
sx q[3];
rz(-1.181447) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(2.7398305) q[2];
sx q[2];
rz(-1.4515452) q[2];
sx q[2];
rz(0.24701992) q[2];
rz(2.5714827) q[3];
sx q[3];
rz(-2.6542122) q[3];
sx q[3];
rz(-2.0232239) q[3];
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
x q[0];
rz(pi/2) q[0];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0982672) q[0];
sx q[0];
rz(-1.7541616) q[0];
sx q[0];
rz(2.7761205) q[0];
rz(0.098516057) q[1];
sx q[1];
rz(-0.64058522) q[1];
sx q[1];
rz(-2.2596901) q[1];
rz(-1.1578887) q[2];
sx q[2];
rz(-1.5831686) q[2];
sx q[2];
rz(-1.8148212) q[2];
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
