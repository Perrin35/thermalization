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
rz(0.61007845) q[0];
sx q[0];
rz(3.371513) q[0];
sx q[0];
rz(11.846677) q[0];
rz(2.5201058) q[1];
sx q[1];
rz(1.4014333) q[1];
sx q[1];
rz(9.2994193) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3925815) q[0];
sx q[0];
rz(-1.1454937) q[0];
sx q[0];
rz(-2.0122737) q[0];
rz(-pi) q[1];
x q[1];
rz(1.4761476) q[2];
sx q[2];
rz(-1.2834335) q[2];
sx q[2];
rz(-0.3269302) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-2.5125864) q[1];
sx q[1];
rz(-3.1401445) q[1];
sx q[1];
rz(1.8667029) q[1];
rz(-pi) q[2];
rz(-2.2940346) q[3];
sx q[3];
rz(-1.8430018) q[3];
sx q[3];
rz(-0.51203007) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(1.3434992) q[2];
sx q[2];
rz(-0.40830475) q[2];
sx q[2];
rz(-2.2957392) q[2];
rz(2.3503303) q[3];
sx q[3];
rz(-0.013412272) q[3];
sx q[3];
rz(-0.066789269) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4132408) q[0];
sx q[0];
rz(-0.46717307) q[0];
sx q[0];
rz(0.043664232) q[0];
rz(1.5751669) q[1];
sx q[1];
rz(-1.3730201) q[1];
sx q[1];
rz(1.498819) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.18066809) q[0];
sx q[0];
rz(-0.60375133) q[0];
sx q[0];
rz(-1.5779499) q[0];
x q[1];
rz(1.5909219) q[2];
sx q[2];
rz(-2.5701036) q[2];
sx q[2];
rz(3.1352941) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(3.1398357) q[1];
sx q[1];
rz(-0.082232177) q[1];
sx q[1];
rz(-2.7539192) q[1];
rz(1.2263857) q[3];
sx q[3];
rz(-1.6020755) q[3];
sx q[3];
rz(2.7117604) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(1.0955536) q[2];
sx q[2];
rz(-2.9914896) q[2];
sx q[2];
rz(-2.6138439) q[2];
rz(2.2857417) q[3];
sx q[3];
rz(-3.1400883) q[3];
sx q[3];
rz(1.9201479) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.39831487) q[0];
sx q[0];
rz(-0.97284955) q[0];
sx q[0];
rz(-1.995218) q[0];
rz(-1.758681) q[1];
sx q[1];
rz(-0.29255602) q[1];
sx q[1];
rz(3.0376099) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2706174) q[0];
sx q[0];
rz(-1.7954602) q[0];
sx q[0];
rz(-1.490834) q[0];
x q[1];
rz(-3.1244794) q[2];
sx q[2];
rz(-1.2984973) q[2];
sx q[2];
rz(0.903331) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.5755363) q[1];
sx q[1];
rz(-0.78278274) q[1];
sx q[1];
rz(2.8498883) q[1];
x q[2];
rz(0.35211925) q[3];
sx q[3];
rz(-2.8330292) q[3];
sx q[3];
rz(-1.7087913) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-2.9581703) q[2];
sx q[2];
rz(-3.1347745) q[2];
sx q[2];
rz(-2.5799694) q[2];
rz(3.0567452) q[3];
sx q[3];
rz(-0.0054797879) q[3];
sx q[3];
rz(-3.1408299) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.7247923) q[0];
sx q[0];
rz(-2.8775207) q[0];
sx q[0];
rz(2.6208139) q[0];
rz(-0.15781038) q[1];
sx q[1];
rz(-0.66693711) q[1];
sx q[1];
rz(-0.075798362) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6822081) q[0];
sx q[0];
rz(-0.59218107) q[0];
sx q[0];
rz(0.20488744) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.5721605) q[2];
sx q[2];
rz(-1.5714386) q[2];
sx q[2];
rz(-0.13880348) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(0.27386777) q[1];
sx q[1];
rz(-2.7368022) q[1];
sx q[1];
rz(-1.936749) q[1];
rz(-pi) q[2];
x q[2];
rz(1.3520207) q[3];
sx q[3];
rz(-0.46275381) q[3];
sx q[3];
rz(0.81172152) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(0.51125222) q[2];
sx q[2];
rz(-3.1255836) q[2];
sx q[2];
rz(-1.7208257) q[2];
rz(-3.1347647) q[3];
sx q[3];
rz(-3.112401) q[3];
sx q[3];
rz(-1.487287) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0066978) q[0];
sx q[0];
rz(-0.20773523) q[0];
sx q[0];
rz(-2.8790706) q[0];
rz(-0.94995704) q[1];
sx q[1];
rz(-0.078153178) q[1];
sx q[1];
rz(0.21771678) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.45633204) q[0];
sx q[0];
rz(-1.963441) q[0];
sx q[0];
rz(3.0079907) q[0];
x q[1];
rz(-1.6545914) q[2];
sx q[2];
rz(-1.4838654) q[2];
sx q[2];
rz(-3.0695335) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(1.3721261) q[1];
sx q[1];
rz(-1.7431152) q[1];
sx q[1];
rz(0.0021688633) q[1];
rz(-pi) q[2];
x q[2];
rz(-3.0721139) q[3];
sx q[3];
rz(-2.0978048) q[3];
sx q[3];
rz(-1.219092) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(0.21506423) q[2];
sx q[2];
rz(-3.1319517) q[2];
sx q[2];
rz(0.24791524) q[2];
rz(1.7667814) q[3];
sx q[3];
rz(-0.050516613) q[3];
sx q[3];
rz(1.4352528) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0404102) q[0];
sx q[0];
rz(-1.8721767) q[0];
sx q[0];
rz(2.5293479) q[0];
rz(2.9712037) q[1];
sx q[1];
rz(-3.0590765) q[1];
sx q[1];
rz(-1.6498529) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9303974) q[0];
sx q[0];
rz(-0.9370417) q[0];
sx q[0];
rz(-0.90715815) q[0];
x q[1];
rz(1.5749919) q[2];
sx q[2];
rz(-1.5853264) q[2];
sx q[2];
rz(-0.57455237) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(2.7266164) q[1];
sx q[1];
rz(-1.6373937) q[1];
sx q[1];
rz(0.29331751) q[1];
rz(-pi) q[2];
rz(-0.18317353) q[3];
sx q[3];
rz(-1.8889931) q[3];
sx q[3];
rz(-0.051200031) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(2.6461569) q[2];
sx q[2];
rz(-0.010224552) q[2];
sx q[2];
rz(-0.23908991) q[2];
rz(2.332989) q[3];
sx q[3];
rz(-0.012152823) q[3];
sx q[3];
rz(-1.3330207) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7072153) q[0];
sx q[0];
rz(-3.0476397) q[0];
sx q[0];
rz(-2.3648426) q[0];
rz(3.1306733) q[1];
sx q[1];
rz(-0.24591406) q[1];
sx q[1];
rz(1.4728665) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4576329) q[0];
sx q[0];
rz(-1.7839971) q[0];
sx q[0];
rz(2.2868566) q[0];
rz(-pi) q[1];
rz(-0.02350925) q[2];
sx q[2];
rz(-1.5644367) q[2];
sx q[2];
rz(-1.3964749) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-0.90172988) q[1];
sx q[1];
rz(-2.3119676) q[1];
sx q[1];
rz(1.4338486) q[1];
rz(-pi) q[2];
rz(0.075447791) q[3];
sx q[3];
rz(-1.8561072) q[3];
sx q[3];
rz(2.4782654) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(1.9933068) q[2];
sx q[2];
rz(-3.1243656) q[2];
sx q[2];
rz(0.41603184) q[2];
rz(2.3038583) q[3];
sx q[3];
rz(-0.13844027) q[3];
sx q[3];
rz(1.4778888) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
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
rz(0.96529043) q[0];
sx q[0];
rz(-0.039881341) q[0];
sx q[0];
rz(0.95529977) q[0];
rz(-1.924986) q[1];
sx q[1];
rz(-2.80426) q[1];
sx q[1];
rz(0.40562707) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.058620166) q[0];
sx q[0];
rz(-0.79917963) q[0];
sx q[0];
rz(2.6891476) q[0];
x q[1];
rz(0.50875278) q[2];
sx q[2];
rz(-0.18157427) q[2];
sx q[2];
rz(-1.5827098) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-0.97598398) q[1];
sx q[1];
rz(-0.17804314) q[1];
sx q[1];
rz(0.44573943) q[1];
x q[2];
rz(-2.9285223) q[3];
sx q[3];
rz(-1.655742) q[3];
sx q[3];
rz(-1.6036878) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(1.8568628) q[2];
sx q[2];
rz(-3.1046107) q[2];
sx q[2];
rz(0.82376087) q[2];
rz(0.13122261) q[3];
sx q[3];
rz(-2.8516912) q[3];
sx q[3];
rz(-2.8141865) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4284978) q[0];
sx q[0];
rz(-0.14398028) q[0];
sx q[0];
rz(2.4289828) q[0];
rz(0.54829848) q[1];
sx q[1];
rz(-0.24990853) q[1];
sx q[1];
rz(-0.20864329) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.945108) q[0];
sx q[0];
rz(-2.5577106) q[0];
sx q[0];
rz(-0.78899011) q[0];
rz(-pi) q[1];
rz(0.018881475) q[2];
sx q[2];
rz(-1.60542) q[2];
sx q[2];
rz(-1.1037301) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-3.123173) q[1];
sx q[1];
rz(-2.5858665) q[1];
sx q[1];
rz(1.5798142) q[1];
rz(-2.3712158) q[3];
sx q[3];
rz(-2.340014) q[3];
sx q[3];
rz(2.1397391) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-1.448552) q[2];
sx q[2];
rz(-0.081144944) q[2];
sx q[2];
rz(2.9244259) q[2];
rz(-2.7740357) q[3];
sx q[3];
rz(-0.03438545) q[3];
sx q[3];
rz(-2.1448081) q[3];
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
sx q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.17598584) q[0];
sx q[0];
rz(-3.0526057) q[0];
sx q[0];
rz(0.17372818) q[0];
rz(1.732775) q[1];
sx q[1];
rz(-1.6830187) q[1];
sx q[1];
rz(-1.6857612) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.834704) q[0];
sx q[0];
rz(-1.406989) q[0];
sx q[0];
rz(0.52229806) q[0];
rz(-pi) q[1];
rz(3.0020724) q[2];
sx q[2];
rz(-1.7762842) q[2];
sx q[2];
rz(0.71256283) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-3.0324754) q[1];
sx q[1];
rz(-0.87504234) q[1];
sx q[1];
rz(-1.4182219) q[1];
x q[2];
rz(-3.0343541) q[3];
sx q[3];
rz(-1.8180366) q[3];
sx q[3];
rz(-3.1013754) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.9659861) q[2];
sx q[2];
rz(-0.49314988) q[2];
sx q[2];
rz(-1.3803587) q[2];
rz(2.4557451) q[3];
sx q[3];
rz(-3.139747) q[3];
sx q[3];
rz(0.67870158) q[3];
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
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3995517) q[0];
sx q[0];
rz(-0.69235943) q[0];
sx q[0];
rz(1.6211744) q[0];
rz(-1.5834658) q[1];
sx q[1];
rz(-1.5065267) q[1];
sx q[1];
rz(-2.9361257) q[1];
rz(-1.4453962) q[2];
sx q[2];
rz(-1.554606) q[2];
sx q[2];
rz(-1.3535624) q[2];
rz(-3.0229983) q[3];
sx q[3];
rz(-1.8684602) q[3];
sx q[3];
rz(2.8475193) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
