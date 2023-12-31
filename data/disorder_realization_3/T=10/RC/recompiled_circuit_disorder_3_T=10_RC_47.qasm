OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(2.2258423) q[0];
sx q[0];
rz(-0.031210829) q[0];
sx q[0];
rz(-0.48506919) q[0];
rz(0.78753161) q[1];
sx q[1];
rz(-1.0163611) q[1];
sx q[1];
rz(2.7273942) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6611377) q[0];
sx q[0];
rz(-1.9248795) q[0];
sx q[0];
rz(-2.7829091) q[0];
rz(1.4235052) q[2];
sx q[2];
rz(-0.5746595) q[2];
sx q[2];
rz(-1.5073843) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(1.011885) q[1];
sx q[1];
rz(-1.9958479) q[1];
sx q[1];
rz(3.0711864) q[1];
x q[2];
rz(-2.7798153) q[3];
sx q[3];
rz(-0.28218111) q[3];
sx q[3];
rz(-0.85229814) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.2177314) q[2];
sx q[2];
rz(-1.8863181) q[2];
sx q[2];
rz(3.1100173) q[2];
rz(1.8850373) q[3];
sx q[3];
rz(-0.50285181) q[3];
sx q[3];
rz(0.52662915) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3027705) q[0];
sx q[0];
rz(-1.4571723) q[0];
sx q[0];
rz(0.1698499) q[0];
rz(-0.70392144) q[1];
sx q[1];
rz(-2.070065) q[1];
sx q[1];
rz(0.53952113) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.3300433) q[0];
sx q[0];
rz(-1.1216315) q[0];
sx q[0];
rz(-0.051785843) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.1621446) q[2];
sx q[2];
rz(-2.2505629) q[2];
sx q[2];
rz(-2.0603927) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.0454228) q[1];
sx q[1];
rz(-1.759937) q[1];
sx q[1];
rz(-2.0778836) q[1];
x q[2];
rz(0.86124729) q[3];
sx q[3];
rz(-1.1303139) q[3];
sx q[3];
rz(1.7064106) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-2.4743621) q[2];
sx q[2];
rz(-1.903406) q[2];
sx q[2];
rz(0.24307069) q[2];
rz(-2.4754751) q[3];
sx q[3];
rz(-0.56454286) q[3];
sx q[3];
rz(1.8977785) q[3];
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
x q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3617525) q[0];
sx q[0];
rz(-3.0267974) q[0];
sx q[0];
rz(-0.4483805) q[0];
rz(-1.7547296) q[1];
sx q[1];
rz(-1.9877537) q[1];
sx q[1];
rz(0.2562491) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8921593) q[0];
sx q[0];
rz(-1.8349577) q[0];
sx q[0];
rz(2.507472) q[0];
rz(-0.44552866) q[2];
sx q[2];
rz(-1.1330714) q[2];
sx q[2];
rz(2.90403) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-0.33830723) q[1];
sx q[1];
rz(-2.3803664) q[1];
sx q[1];
rz(2.0558946) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.9647397) q[3];
sx q[3];
rz(-1.5336509) q[3];
sx q[3];
rz(-1.1249441) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.8024575) q[2];
sx q[2];
rz(-0.70636237) q[2];
sx q[2];
rz(0.57470542) q[2];
rz(1.3556708) q[3];
sx q[3];
rz(-1.971743) q[3];
sx q[3];
rz(-1.2333966) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.213585) q[0];
sx q[0];
rz(-1.428823) q[0];
sx q[0];
rz(2.8821049) q[0];
rz(-1.150594) q[1];
sx q[1];
rz(-1.3535627) q[1];
sx q[1];
rz(-0.73192275) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.61011945) q[0];
sx q[0];
rz(-1.1024794) q[0];
sx q[0];
rz(2.5584695) q[0];
rz(0.133693) q[2];
sx q[2];
rz(-0.72172726) q[2];
sx q[2];
rz(2.0267817) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-0.54938984) q[1];
sx q[1];
rz(-1.376774) q[1];
sx q[1];
rz(-2.8787896) q[1];
x q[2];
rz(-2.9210806) q[3];
sx q[3];
rz(-1.6933428) q[3];
sx q[3];
rz(1.3361271) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(1.6966454) q[2];
sx q[2];
rz(-1.8680633) q[2];
sx q[2];
rz(2.990492) q[2];
rz(0.54667306) q[3];
sx q[3];
rz(-1.0497382) q[3];
sx q[3];
rz(0.37724075) q[3];
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
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5916409) q[0];
sx q[0];
rz(-2.0786091) q[0];
sx q[0];
rz(1.2623825) q[0];
rz(1.6732015) q[1];
sx q[1];
rz(-0.60931283) q[1];
sx q[1];
rz(-0.79777065) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7724458) q[0];
sx q[0];
rz(-1.6909084) q[0];
sx q[0];
rz(3.1390879) q[0];
rz(-pi) q[1];
rz(-1.7902137) q[2];
sx q[2];
rz(-1.2741538) q[2];
sx q[2];
rz(2.9079633) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-3.0610173) q[1];
sx q[1];
rz(-1.4772381) q[1];
sx q[1];
rz(2.1993125) q[1];
x q[2];
rz(-1.9005152) q[3];
sx q[3];
rz(-1.1493249) q[3];
sx q[3];
rz(0.7082522) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.795934) q[2];
sx q[2];
rz(-0.63085932) q[2];
sx q[2];
rz(2.8395555) q[2];
rz(1.1473514) q[3];
sx q[3];
rz(-1.6796422) q[3];
sx q[3];
rz(2.5938477) q[3];
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
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.40925947) q[0];
sx q[0];
rz(-3.0497666) q[0];
sx q[0];
rz(-1.1557895) q[0];
rz(-2.0571158) q[1];
sx q[1];
rz(-2.1613354) q[1];
sx q[1];
rz(-3.0715122) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.72567155) q[0];
sx q[0];
rz(-1.7416746) q[0];
sx q[0];
rz(0.15539774) q[0];
rz(-pi) q[1];
rz(-2.1224535) q[2];
sx q[2];
rz(-0.86691228) q[2];
sx q[2];
rz(0.64955074) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-1.8225704) q[1];
sx q[1];
rz(-2.550107) q[1];
sx q[1];
rz(2.849008) q[1];
rz(-0.035404215) q[3];
sx q[3];
rz(-1.6453214) q[3];
sx q[3];
rz(-2.7217334) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-2.8391116) q[2];
sx q[2];
rz(-1.182686) q[2];
sx q[2];
rz(-0.62136674) q[2];
rz(-1.4403884) q[3];
sx q[3];
rz(-2.6337603) q[3];
sx q[3];
rz(2.5682209) q[3];
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
sx q[0];
rz(-pi) q[0];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0550585) q[0];
sx q[0];
rz(-1.4165514) q[0];
sx q[0];
rz(-2.281718) q[0];
rz(-1.2043918) q[1];
sx q[1];
rz(-2.269373) q[1];
sx q[1];
rz(-3.133657) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7746349) q[0];
sx q[0];
rz(-1.3152221) q[0];
sx q[0];
rz(-0.5704244) q[0];
rz(-pi) q[1];
x q[1];
rz(2.7107312) q[2];
sx q[2];
rz(-2.2164946) q[2];
sx q[2];
rz(2.8889887) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-1.7759526) q[1];
sx q[1];
rz(-1.0675634) q[1];
sx q[1];
rz(0.788049) q[1];
x q[2];
rz(2.2357335) q[3];
sx q[3];
rz(-2.6568036) q[3];
sx q[3];
rz(-1.0672027) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-0.69592151) q[2];
sx q[2];
rz(-1.9182703) q[2];
sx q[2];
rz(0.41637862) q[2];
rz(-1.773206) q[3];
sx q[3];
rz(-1.2973283) q[3];
sx q[3];
rz(-0.90464512) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
x q[3];
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
rz(0.96173441) q[0];
sx q[0];
rz(-2.8680153) q[0];
sx q[0];
rz(-2.7767048) q[0];
rz(-2.2015613) q[1];
sx q[1];
rz(-2.6049728) q[1];
sx q[1];
rz(1.6392802) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.72909268) q[0];
sx q[0];
rz(-1.3206498) q[0];
sx q[0];
rz(-1.5604707) q[0];
x q[1];
rz(2.084923) q[2];
sx q[2];
rz(-0.80386111) q[2];
sx q[2];
rz(0.67509292) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(2.9598436) q[1];
sx q[1];
rz(-1.0105003) q[1];
sx q[1];
rz(0.77397857) q[1];
rz(-pi) q[2];
rz(-1.7765462) q[3];
sx q[3];
rz(-1.660941) q[3];
sx q[3];
rz(0.84116018) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(0.49729785) q[2];
sx q[2];
rz(-2.6361894) q[2];
sx q[2];
rz(-2.0765182) q[2];
rz(-2.8403357) q[3];
sx q[3];
rz(-1.732429) q[3];
sx q[3];
rz(-1.6206954) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.97312462) q[0];
sx q[0];
rz(-0.081806101) q[0];
sx q[0];
rz(2.705943) q[0];
rz(-1.3849974) q[1];
sx q[1];
rz(-2.6711617) q[1];
sx q[1];
rz(0.41697821) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1736261) q[0];
sx q[0];
rz(-2.2244503) q[0];
sx q[0];
rz(-3.0941512) q[0];
rz(-pi) q[1];
x q[1];
rz(0.99879361) q[2];
sx q[2];
rz(-0.98888328) q[2];
sx q[2];
rz(-1.7172161) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-1.4638085) q[1];
sx q[1];
rz(-0.52075547) q[1];
sx q[1];
rz(-1.2893454) q[1];
rz(-pi) q[2];
x q[2];
rz(2.3973893) q[3];
sx q[3];
rz(-0.19072285) q[3];
sx q[3];
rz(2.423559) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-0.10432648) q[2];
sx q[2];
rz(-1.4929079) q[2];
sx q[2];
rz(2.9157675) q[2];
rz(-2.9337692) q[3];
sx q[3];
rz(-2.4184629) q[3];
sx q[3];
rz(-0.58661714) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
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
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.41480961) q[0];
sx q[0];
rz(-2.8746958) q[0];
sx q[0];
rz(1.6171932) q[0];
rz(-2.1879451) q[1];
sx q[1];
rz(-1.8667659) q[1];
sx q[1];
rz(-1.8189925) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.029862558) q[0];
sx q[0];
rz(-2.7435281) q[0];
sx q[0];
rz(1.8882621) q[0];
rz(-1.1216713) q[2];
sx q[2];
rz(-1.7977062) q[2];
sx q[2];
rz(-0.34044468) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(3.0364914) q[1];
sx q[1];
rz(-2.0836012) q[1];
sx q[1];
rz(-0.75863104) q[1];
rz(1.7694468) q[3];
sx q[3];
rz(-2.47654) q[3];
sx q[3];
rz(-2.8952451) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-1.7913197) q[2];
sx q[2];
rz(-2.0952756) q[2];
sx q[2];
rz(2.1255169) q[2];
rz(-1.2223876) q[3];
sx q[3];
rz(-0.17761579) q[3];
sx q[3];
rz(0.55541742) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
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
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9983457) q[0];
sx q[0];
rz(-2.2194942) q[0];
sx q[0];
rz(1.0794328) q[0];
rz(-1.7779508) q[1];
sx q[1];
rz(-1.9121871) q[1];
sx q[1];
rz(-1.8180064) q[1];
rz(-3.1238363) q[2];
sx q[2];
rz(-2.6323071) q[2];
sx q[2];
rz(1.4321362) q[2];
rz(1.2033403) q[3];
sx q[3];
rz(-0.24693476) q[3];
sx q[3];
rz(0.91528391) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
