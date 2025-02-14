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
rz(4.3920565) q[1];
sx q[1];
rz(5.1345706) q[1];
sx q[1];
rz(13.486025) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5772668) q[0];
sx q[0];
rz(-2.3503605) q[0];
sx q[0];
rz(-1.7971695) q[0];
rz(-pi) q[1];
rz(0.55962015) q[2];
sx q[2];
rz(-2.3336448) q[2];
sx q[2];
rz(1.9576278) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.8383467) q[1];
sx q[1];
rz(-0.64323264) q[1];
sx q[1];
rz(0.077416181) q[1];
rz(-pi) q[2];
rz(-1.7149431) q[3];
sx q[3];
rz(-2.1021772) q[3];
sx q[3];
rz(1.8559769) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-0.93996843) q[2];
sx q[2];
rz(-2.2087966) q[2];
sx q[2];
rz(-1.0533818) q[2];
rz(0.71422226) q[3];
sx q[3];
rz(-0.37319365) q[3];
sx q[3];
rz(-1.8478954) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6222222) q[0];
sx q[0];
rz(-0.70835963) q[0];
sx q[0];
rz(-0.57193065) q[0];
rz(1.8427303) q[1];
sx q[1];
rz(-0.79890257) q[1];
sx q[1];
rz(2.3518708) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.1173218) q[0];
sx q[0];
rz(-2.824475) q[0];
sx q[0];
rz(-1.8834524) q[0];
rz(-pi) q[1];
x q[1];
rz(2.9339173) q[2];
sx q[2];
rz(-2.4768314) q[2];
sx q[2];
rz(1.1663811) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(1.4724628) q[1];
sx q[1];
rz(-1.8505368) q[1];
sx q[1];
rz(-0.6196687) q[1];
rz(-pi) q[2];
rz(-2.4667707) q[3];
sx q[3];
rz(-1.9673507) q[3];
sx q[3];
rz(0.048585437) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(2.2907437) q[2];
sx q[2];
rz(-0.76883832) q[2];
sx q[2];
rz(1.7791465) q[2];
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
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0027851) q[0];
sx q[0];
rz(-1.0951575) q[0];
sx q[0];
rz(-2.441067) q[0];
rz(-0.48577148) q[1];
sx q[1];
rz(-2.3120717) q[1];
sx q[1];
rz(2.9170759) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1248847) q[0];
sx q[0];
rz(-1.9210235) q[0];
sx q[0];
rz(1.7279529) q[0];
x q[1];
rz(3.0430138) q[2];
sx q[2];
rz(-2.1469627) q[2];
sx q[2];
rz(-2.2599932) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-0.91229328) q[1];
sx q[1];
rz(-2.5535431) q[1];
sx q[1];
rz(-1.4347784) q[1];
rz(1.6979917) q[3];
sx q[3];
rz(-2.4089097) q[3];
sx q[3];
rz(-0.81058433) q[3];
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
rz(-1.1578355) q[3];
sx q[3];
rz(2.4750347) q[3];
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
rz(-pi) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9662358) q[0];
sx q[0];
rz(-0.55713621) q[0];
sx q[0];
rz(0.019158451) q[0];
rz(1.3592023) q[1];
sx q[1];
rz(-1.7117585) q[1];
sx q[1];
rz(-3.137099) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2990103) q[0];
sx q[0];
rz(-2.0662464) q[0];
sx q[0];
rz(-0.61758496) q[0];
x q[1];
rz(-2.0199297) q[2];
sx q[2];
rz(-2.0259078) q[2];
sx q[2];
rz(-0.016985026) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(3.0862866) q[1];
sx q[1];
rz(-1.9290566) q[1];
sx q[1];
rz(-2.9614224) q[1];
rz(-pi) q[2];
rz(2.717797) q[3];
sx q[3];
rz(-2.9188699) q[3];
sx q[3];
rz(0.84171695) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-0.4517453) q[2];
sx q[2];
rz(-1.0901901) q[2];
sx q[2];
rz(1.4534265) q[2];
rz(-1.8870185) q[3];
sx q[3];
rz(-1.5213608) q[3];
sx q[3];
rz(-2.5793251) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3314064) q[0];
sx q[0];
rz(-1.2368546) q[0];
sx q[0];
rz(1.1619262) q[0];
rz(2.3792073) q[1];
sx q[1];
rz(-0.98680174) q[1];
sx q[1];
rz(0.42957482) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7976101) q[0];
sx q[0];
rz(-1.5344041) q[0];
sx q[0];
rz(1.9790566) q[0];
rz(-pi) q[1];
x q[1];
rz(-7/(8*pi)) q[2];
sx q[2];
rz(-2.3633011) q[2];
sx q[2];
rz(-2.8964004) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-2.4641494) q[1];
sx q[1];
rz(-0.92532571) q[1];
sx q[1];
rz(0.42663017) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.91752538) q[3];
sx q[3];
rz(-1.297373) q[3];
sx q[3];
rz(-1.3454252) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(0.21542159) q[2];
sx q[2];
rz(-1.8502219) q[2];
sx q[2];
rz(1.8298836) q[2];
rz(2.6808776) q[3];
sx q[3];
rz(-1.0034794) q[3];
sx q[3];
rz(-3.0084012) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8120414) q[0];
sx q[0];
rz(-2.9819745) q[0];
sx q[0];
rz(1.106369) q[0];
rz(-3.0270992) q[1];
sx q[1];
rz(-1.0527) q[1];
sx q[1];
rz(3.0879424) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.69218) q[0];
sx q[0];
rz(-0.25568889) q[0];
sx q[0];
rz(-1.5480255) q[0];
x q[1];
rz(-1.6752376) q[2];
sx q[2];
rz(-0.76600961) q[2];
sx q[2];
rz(-3.0867689) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(2.1785894) q[1];
sx q[1];
rz(-0.51437639) q[1];
sx q[1];
rz(1.4609171) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.63416449) q[3];
sx q[3];
rz(-2.0863219) q[3];
sx q[3];
rz(-0.17532119) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.8657118) q[2];
sx q[2];
rz(-1.2578473) q[2];
sx q[2];
rz(-0.50714058) q[2];
rz(1.7670613) q[3];
sx q[3];
rz(-1.5406698) q[3];
sx q[3];
rz(-1.1486294) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5116665) q[0];
sx q[0];
rz(-2.0138795) q[0];
sx q[0];
rz(3.1003057) q[0];
rz(-0.73829007) q[1];
sx q[1];
rz(-1.2246882) q[1];
sx q[1];
rz(-0.98947492) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8007198) q[0];
sx q[0];
rz(-1.1015019) q[0];
sx q[0];
rz(-1.0020026) q[0];
rz(2.6648952) q[2];
sx q[2];
rz(-0.39604353) q[2];
sx q[2];
rz(-0.11644289) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-1.3937903) q[1];
sx q[1];
rz(-0.89266864) q[1];
sx q[1];
rz(-0.80121471) q[1];
rz(-pi) q[2];
rz(-2.2253898) q[3];
sx q[3];
rz(-1.8713017) q[3];
sx q[3];
rz(0.88778561) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.5993293) q[2];
sx q[2];
rz(-1.0487391) q[2];
sx q[2];
rz(-0.85912022) q[2];
rz(0.011693444) q[3];
sx q[3];
rz(-1.499736) q[3];
sx q[3];
rz(0.11277994) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.82666731) q[0];
sx q[0];
rz(-2.5496917) q[0];
sx q[0];
rz(2.9097606) q[0];
rz(-2.5595317) q[1];
sx q[1];
rz(-1.3287909) q[1];
sx q[1];
rz(1.9225165) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.40755475) q[0];
sx q[0];
rz(-2.1269607) q[0];
sx q[0];
rz(-1.708605) q[0];
rz(-pi) q[1];
x q[1];
rz(2.788576) q[2];
sx q[2];
rz(-1.5215877) q[2];
sx q[2];
rz(-0.69248688) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-2.4640567) q[1];
sx q[1];
rz(-1.4194064) q[1];
sx q[1];
rz(1.7580126) q[1];
x q[2];
rz(-0.7557423) q[3];
sx q[3];
rz(-1.983194) q[3];
sx q[3];
rz(2.2568767) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(1.9570738) q[2];
sx q[2];
rz(-1.7902713) q[2];
sx q[2];
rz(-0.78682023) q[2];
rz(-1.5015073) q[3];
sx q[3];
rz(-0.4314751) q[3];
sx q[3];
rz(-1.7155581) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.96104923) q[0];
sx q[0];
rz(-1.7683872) q[0];
sx q[0];
rz(-2.2013262) q[0];
rz(-0.44479784) q[1];
sx q[1];
rz(-0.85591379) q[1];
sx q[1];
rz(0.14642265) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.1347552) q[0];
sx q[0];
rz(-1.2723347) q[0];
sx q[0];
rz(1.439007) q[0];
rz(-pi) q[1];
rz(1.8249885) q[2];
sx q[2];
rz(-1.3066402) q[2];
sx q[2];
rz(0.025042695) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(2.2487029) q[1];
sx q[1];
rz(-1.8161402) q[1];
sx q[1];
rz(-0.69458436) q[1];
rz(-0.10036631) q[3];
sx q[3];
rz(-0.79381889) q[3];
sx q[3];
rz(0.19089261) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(0.08708295) q[2];
sx q[2];
rz(-1.4887709) q[2];
sx q[2];
rz(0.83756891) q[2];
rz(-2.2086823) q[3];
sx q[3];
rz(-0.98740238) q[3];
sx q[3];
rz(2.3605409) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
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
rz(-2.8785716) q[0];
sx q[0];
rz(-1.2282547) q[0];
sx q[0];
rz(1.3680869) q[0];
rz(0.076016501) q[1];
sx q[1];
rz(-0.58034211) q[1];
sx q[1];
rz(-2.0215633) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.66066) q[0];
sx q[0];
rz(-1.8802065) q[0];
sx q[0];
rz(-2.8641228) q[0];
rz(-pi) q[1];
x q[1];
rz(1.8755959) q[2];
sx q[2];
rz(-2.2561361) q[2];
sx q[2];
rz(0.03396578) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-3.0894153) q[1];
sx q[1];
rz(-1.4751768) q[1];
sx q[1];
rz(-1.862942) q[1];
rz(3.038123) q[3];
sx q[3];
rz(-1.4082239) q[3];
sx q[3];
rz(1.9601456) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(0.40176216) q[2];
sx q[2];
rz(-1.4515452) q[2];
sx q[2];
rz(-0.24701992) q[2];
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
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0433255) q[0];
sx q[0];
rz(-1.7541616) q[0];
sx q[0];
rz(2.7761205) q[0];
rz(-0.098516057) q[1];
sx q[1];
rz(-2.5010074) q[1];
sx q[1];
rz(0.88190257) q[1];
rz(-0.013507387) q[2];
sx q[2];
rz(-1.9836704) q[2];
sx q[2];
rz(-0.23860485) q[2];
rz(1.3104812) q[3];
sx q[3];
rz(-1.9394939) q[3];
sx q[3];
rz(0.60422411) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
