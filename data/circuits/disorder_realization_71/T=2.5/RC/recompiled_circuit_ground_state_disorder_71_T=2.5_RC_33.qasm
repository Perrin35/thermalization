OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-2.3040721) q[0];
sx q[0];
rz(-2.4763595) q[0];
sx q[0];
rz(-1.2256149) q[0];
rz(2.4576814) q[1];
sx q[1];
rz(-0.39165762) q[1];
sx q[1];
rz(-2.439523) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.41034979) q[0];
sx q[0];
rz(-0.87603891) q[0];
sx q[0];
rz(1.4175416) q[0];
rz(-pi) q[1];
x q[1];
rz(1.6984387) q[2];
sx q[2];
rz(-0.1662456) q[2];
sx q[2];
rz(0.12285168) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-2.7154123) q[1];
sx q[1];
rz(-0.75580929) q[1];
sx q[1];
rz(-2.4017357) q[1];
rz(-2.2625173) q[3];
sx q[3];
rz(-1.4180776) q[3];
sx q[3];
rz(-0.16593753) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.94032732) q[2];
sx q[2];
rz(-1.5387115) q[2];
sx q[2];
rz(-2.8931457) q[2];
rz(-2.0072319) q[3];
sx q[3];
rz(-3.0076707) q[3];
sx q[3];
rz(1.1321446) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.092875384) q[0];
sx q[0];
rz(-2.5889914) q[0];
sx q[0];
rz(-1.5930814) q[0];
rz(-2.6379207) q[1];
sx q[1];
rz(-1.2232989) q[1];
sx q[1];
rz(-0.37685397) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5232613) q[0];
sx q[0];
rz(-1.0987765) q[0];
sx q[0];
rz(1.5923772) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.3289248) q[2];
sx q[2];
rz(-2.118131) q[2];
sx q[2];
rz(-0.081204435) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.845692) q[1];
sx q[1];
rz(-1.0949425) q[1];
sx q[1];
rz(-1.8527387) q[1];
x q[2];
rz(0.33263388) q[3];
sx q[3];
rz(-1.4983699) q[3];
sx q[3];
rz(0.70055947) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-2.5571931) q[2];
sx q[2];
rz(-0.69584766) q[2];
sx q[2];
rz(2.2759571) q[2];
rz(1.4731167) q[3];
sx q[3];
rz(-2.1560463) q[3];
sx q[3];
rz(-2.1540811) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5675548) q[0];
sx q[0];
rz(-1.8542629) q[0];
sx q[0];
rz(1.0512742) q[0];
rz(1.4569262) q[1];
sx q[1];
rz(-1.8345865) q[1];
sx q[1];
rz(-1.5164703) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.595436) q[0];
sx q[0];
rz(-1.8519028) q[0];
sx q[0];
rz(-0.77451046) q[0];
rz(-pi) q[1];
x q[1];
rz(2.3651354) q[2];
sx q[2];
rz(-2.2946022) q[2];
sx q[2];
rz(3.0645097) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(0.96788266) q[1];
sx q[1];
rz(-2.5213402) q[1];
sx q[1];
rz(-2.2642676) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.72318913) q[3];
sx q[3];
rz(-0.08130493) q[3];
sx q[3];
rz(-1.7618084) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.7463344) q[2];
sx q[2];
rz(-1.0552152) q[2];
sx q[2];
rz(-0.22221097) q[2];
rz(2.639751) q[3];
sx q[3];
rz(-1.7343438) q[3];
sx q[3];
rz(2.3327904) q[3];
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
x q[0];
rz(-pi/2) q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2489081) q[0];
sx q[0];
rz(-2.6539256) q[0];
sx q[0];
rz(1.9768313) q[0];
rz(-2.248863) q[1];
sx q[1];
rz(-1.3803394) q[1];
sx q[1];
rz(-2.8597615) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4455216) q[0];
sx q[0];
rz(-2.0508678) q[0];
sx q[0];
rz(1.4351373) q[0];
rz(-pi) q[1];
rz(-2.495419) q[2];
sx q[2];
rz(-0.82662383) q[2];
sx q[2];
rz(-1.7648362) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-2.6161595) q[1];
sx q[1];
rz(-1.6405917) q[1];
sx q[1];
rz(-1.4074123) q[1];
rz(-pi) q[2];
rz(-1.6244333) q[3];
sx q[3];
rz(-2.6948018) q[3];
sx q[3];
rz(0.051366816) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(2.0622327) q[2];
sx q[2];
rz(-0.56228176) q[2];
sx q[2];
rz(-2.6083561) q[2];
rz(-1.0821292) q[3];
sx q[3];
rz(-2.5366668) q[3];
sx q[3];
rz(2.3281039) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.37394062) q[0];
sx q[0];
rz(-0.39893183) q[0];
sx q[0];
rz(1.3237413) q[0];
rz(2.3228877) q[1];
sx q[1];
rz(-2.3981514) q[1];
sx q[1];
rz(2.0735819) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0518347) q[0];
sx q[0];
rz(-1.7300419) q[0];
sx q[0];
rz(1.6510886) q[0];
x q[1];
rz(2.2862412) q[2];
sx q[2];
rz(-0.87424874) q[2];
sx q[2];
rz(-1.0698505) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(0.45324486) q[1];
sx q[1];
rz(-1.6736341) q[1];
sx q[1];
rz(2.4842841) q[1];
rz(2.976494) q[3];
sx q[3];
rz(-2.0378651) q[3];
sx q[3];
rz(1.1822753) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(1.9615122) q[2];
sx q[2];
rz(-2.027812) q[2];
sx q[2];
rz(0.86165825) q[2];
rz(-2.1152451) q[3];
sx q[3];
rz(-1.9828826) q[3];
sx q[3];
rz(-2.0238743) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3468129) q[0];
sx q[0];
rz(-0.81951278) q[0];
sx q[0];
rz(-2.2487707) q[0];
rz(-0.18094856) q[1];
sx q[1];
rz(-2.4632958) q[1];
sx q[1];
rz(-0.13042626) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.82662) q[0];
sx q[0];
rz(-0.19222799) q[0];
sx q[0];
rz(-2.1703224) q[0];
rz(-pi) q[1];
rz(2.0753292) q[2];
sx q[2];
rz(-1.0257402) q[2];
sx q[2];
rz(1.2912599) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(2.0196258) q[1];
sx q[1];
rz(-1.5415915) q[1];
sx q[1];
rz(-1.8080416) q[1];
x q[2];
rz(2.942286) q[3];
sx q[3];
rz(-2.350507) q[3];
sx q[3];
rz(1.6419322) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-0.22455198) q[2];
sx q[2];
rz(-1.9150534) q[2];
sx q[2];
rz(-0.87635931) q[2];
rz(-1.2419491) q[3];
sx q[3];
rz(-0.94049898) q[3];
sx q[3];
rz(2.131264) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7414311) q[0];
sx q[0];
rz(-3.042996) q[0];
sx q[0];
rz(2.7778991) q[0];
rz(2.3997276) q[1];
sx q[1];
rz(-1.0262841) q[1];
sx q[1];
rz(2.738764) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9972853) q[0];
sx q[0];
rz(-1.3308405) q[0];
sx q[0];
rz(-1.994094) q[0];
x q[1];
rz(1.9697519) q[2];
sx q[2];
rz(-1.9213994) q[2];
sx q[2];
rz(0.95870852) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(1.6322569) q[1];
sx q[1];
rz(-1.4493353) q[1];
sx q[1];
rz(3.0640825) q[1];
rz(-pi) q[2];
x q[2];
rz(1.8419344) q[3];
sx q[3];
rz(-0.95356546) q[3];
sx q[3];
rz(0.24880508) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(1.3591298) q[2];
sx q[2];
rz(-2.7733347) q[2];
sx q[2];
rz(-1.5472319) q[2];
rz(-0.83526978) q[3];
sx q[3];
rz(-2.1850977) q[3];
sx q[3];
rz(2.6529151) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.35992026) q[0];
sx q[0];
rz(-1.3363573) q[0];
sx q[0];
rz(-2.6265327) q[0];
rz(-1.7022279) q[1];
sx q[1];
rz(-1.2934338) q[1];
sx q[1];
rz(0.39915592) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3278807) q[0];
sx q[0];
rz(-1.7017168) q[0];
sx q[0];
rz(1.5689724) q[0];
rz(-pi) q[1];
rz(-2.299526) q[2];
sx q[2];
rz(-1.3883603) q[2];
sx q[2];
rz(1.6853756) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(1.5134352) q[1];
sx q[1];
rz(-1.6353459) q[1];
sx q[1];
rz(1.8316395) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.5695482) q[3];
sx q[3];
rz(-1.5051418) q[3];
sx q[3];
rz(-1.9492689) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(1.9814375) q[2];
sx q[2];
rz(-1.258054) q[2];
sx q[2];
rz(-1.0464279) q[2];
rz(0.6662755) q[3];
sx q[3];
rz(-2.8847238) q[3];
sx q[3];
rz(1.0506786) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1995354) q[0];
sx q[0];
rz(-0.10534795) q[0];
sx q[0];
rz(0.60607213) q[0];
rz(2.3145158) q[1];
sx q[1];
rz(-1.8111633) q[1];
sx q[1];
rz(-1.3333295) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5510941) q[0];
sx q[0];
rz(-2.7685389) q[0];
sx q[0];
rz(0.28753745) q[0];
rz(-pi) q[1];
x q[1];
rz(2.016388) q[2];
sx q[2];
rz(-2.2632709) q[2];
sx q[2];
rz(-2.8852579) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(0.096790678) q[1];
sx q[1];
rz(-1.6875182) q[1];
sx q[1];
rz(0.90116426) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.656084) q[3];
sx q[3];
rz(-1.1138289) q[3];
sx q[3];
rz(2.5449139) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(2.1313608) q[2];
sx q[2];
rz(-2.2183245) q[2];
sx q[2];
rz(-2.7566747) q[2];
rz(-2.5108003) q[3];
sx q[3];
rz(-1.8600978) q[3];
sx q[3];
rz(1.3425286) q[3];
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
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.25403062) q[0];
sx q[0];
rz(-1.7429202) q[0];
sx q[0];
rz(-1.4400462) q[0];
rz(2.4002659) q[1];
sx q[1];
rz(-1.7404375) q[1];
sx q[1];
rz(-2.8828566) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.16782204) q[0];
sx q[0];
rz(-0.24594618) q[0];
sx q[0];
rz(1.587338) q[0];
rz(1.0503631) q[2];
sx q[2];
rz(-2.5602617) q[2];
sx q[2];
rz(2.5985825) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-2.179271) q[1];
sx q[1];
rz(-1.2291161) q[1];
sx q[1];
rz(0.82247295) q[1];
x q[2];
rz(-2.7388938) q[3];
sx q[3];
rz(-0.81953632) q[3];
sx q[3];
rz(2.0961026) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.4197293) q[2];
sx q[2];
rz(-1.1121007) q[2];
sx q[2];
rz(-1.3191684) q[2];
rz(-2.3223274) q[3];
sx q[3];
rz(-2.0115439) q[3];
sx q[3];
rz(-0.54660249) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7405613) q[0];
sx q[0];
rz(-0.87420976) q[0];
sx q[0];
rz(2.5921205) q[0];
rz(1.7026547) q[1];
sx q[1];
rz(-2.4733652) q[1];
sx q[1];
rz(0.53818902) q[1];
rz(-0.010942608) q[2];
sx q[2];
rz(-0.86075114) q[2];
sx q[2];
rz(1.1759947) q[2];
rz(1.9520252) q[3];
sx q[3];
rz(-0.59567957) q[3];
sx q[3];
rz(2.4761562) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
