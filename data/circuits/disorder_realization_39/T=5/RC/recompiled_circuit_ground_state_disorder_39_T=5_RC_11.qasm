OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(0.47867632) q[0];
sx q[0];
rz(7.8362099) q[0];
sx q[0];
rz(10.683164) q[0];
rz(1.2619184) q[1];
sx q[1];
rz(-2.6231397) q[1];
sx q[1];
rz(-0.060103091) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3451884) q[0];
sx q[0];
rz(-3.0462) q[0];
sx q[0];
rz(-0.32890396) q[0];
rz(-pi) q[1];
rz(-1.41729) q[2];
sx q[2];
rz(-2.058284) q[2];
sx q[2];
rz(1.1660898) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.2138252) q[1];
sx q[1];
rz(-2.7672184) q[1];
sx q[1];
rz(-1.8273749) q[1];
rz(-pi) q[2];
x q[2];
rz(2.6782436) q[3];
sx q[3];
rz(-1.0245205) q[3];
sx q[3];
rz(2.0093371) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-0.95187676) q[2];
sx q[2];
rz(-2.3592301) q[2];
sx q[2];
rz(-0.18048364) q[2];
rz(2.8372724) q[3];
sx q[3];
rz(-0.92420998) q[3];
sx q[3];
rz(2.9144104) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.75905281) q[0];
sx q[0];
rz(-2.5681684) q[0];
sx q[0];
rz(-3.0376814) q[0];
rz(-0.78481627) q[1];
sx q[1];
rz(-1.8321313) q[1];
sx q[1];
rz(-2.5596502) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.7933998) q[0];
sx q[0];
rz(-1.5485067) q[0];
sx q[0];
rz(-0.53198708) q[0];
rz(-2.632276) q[2];
sx q[2];
rz(-1.881372) q[2];
sx q[2];
rz(-2.3485225) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-0.84060639) q[1];
sx q[1];
rz(-2.0475629) q[1];
sx q[1];
rz(-0.029164) q[1];
x q[2];
rz(-2.7304303) q[3];
sx q[3];
rz(-1.9389279) q[3];
sx q[3];
rz(1.0228093) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-2.4265784) q[2];
sx q[2];
rz(-2.7535186) q[2];
sx q[2];
rz(-1.0283872) q[2];
rz(-2.5905124) q[3];
sx q[3];
rz(-1.9469399) q[3];
sx q[3];
rz(-0.50857956) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
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
rz(-0.55117637) q[0];
sx q[0];
rz(-1.6352147) q[0];
sx q[0];
rz(-1.8835541) q[0];
rz(-2.538077) q[1];
sx q[1];
rz(-1.8305093) q[1];
sx q[1];
rz(2.9579732) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.98417898) q[0];
sx q[0];
rz(-1.8350198) q[0];
sx q[0];
rz(1.9963032) q[0];
rz(-pi) q[1];
x q[1];
rz(1.7043978) q[2];
sx q[2];
rz(-1.0935942) q[2];
sx q[2];
rz(1.5988775) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(1.2752876) q[1];
sx q[1];
rz(-0.86967378) q[1];
sx q[1];
rz(2.5661468) q[1];
rz(2.9345064) q[3];
sx q[3];
rz(-1.5303751) q[3];
sx q[3];
rz(-0.8881027) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-0.52643481) q[2];
sx q[2];
rz(-1.7650471) q[2];
sx q[2];
rz(-0.60205013) q[2];
rz(2.708882) q[3];
sx q[3];
rz(-1.5940462) q[3];
sx q[3];
rz(1.6656779) q[3];
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
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3848307) q[0];
sx q[0];
rz(-2.6669406) q[0];
sx q[0];
rz(-2.5153644) q[0];
rz(1.2681883) q[1];
sx q[1];
rz(-1.7363997) q[1];
sx q[1];
rz(2.7558806) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.69067467) q[0];
sx q[0];
rz(-1.1407778) q[0];
sx q[0];
rz(-2.3818092) q[0];
rz(-pi) q[1];
x q[1];
rz(0.63718225) q[2];
sx q[2];
rz(-2.4546625) q[2];
sx q[2];
rz(-2.1241344) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(1.6059162) q[1];
sx q[1];
rz(-1.7711346) q[1];
sx q[1];
rz(-1.5027593) q[1];
x q[2];
rz(-3.0195531) q[3];
sx q[3];
rz(-1.0309848) q[3];
sx q[3];
rz(0.5468747) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-1.4422153) q[2];
sx q[2];
rz(-1.135004) q[2];
sx q[2];
rz(-2.8507612) q[2];
rz(-1.95131) q[3];
sx q[3];
rz(-0.61605993) q[3];
sx q[3];
rz(-2.0975838) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6749343) q[0];
sx q[0];
rz(-3.0789154) q[0];
sx q[0];
rz(-1.4516996) q[0];
rz(2.2845204) q[1];
sx q[1];
rz(-1.2985726) q[1];
sx q[1];
rz(0.42795408) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.33855406) q[0];
sx q[0];
rz(-2.0074275) q[0];
sx q[0];
rz(3.0019041) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.5922115) q[2];
sx q[2];
rz(-0.81839389) q[2];
sx q[2];
rz(2.9894418) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-2.1623858) q[1];
sx q[1];
rz(-0.89079327) q[1];
sx q[1];
rz(2.0335781) q[1];
rz(-pi) q[2];
rz(2.487872) q[3];
sx q[3];
rz(-0.80721426) q[3];
sx q[3];
rz(-1.2098055) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(1.2032623) q[2];
sx q[2];
rz(-2.0945956) q[2];
sx q[2];
rz(-2.9948575) q[2];
rz(0.19715582) q[3];
sx q[3];
rz(-1.6517755) q[3];
sx q[3];
rz(2.3728235) q[3];
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
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.41106221) q[0];
sx q[0];
rz(-1.0163607) q[0];
sx q[0];
rz(-2.9732669) q[0];
rz(0.39237818) q[1];
sx q[1];
rz(-2.7057251) q[1];
sx q[1];
rz(1.4515152) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0453059) q[0];
sx q[0];
rz(-1.1340965) q[0];
sx q[0];
rz(0.023604579) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.8658357) q[2];
sx q[2];
rz(-2.3988798) q[2];
sx q[2];
rz(-0.64464049) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(2.299356) q[1];
sx q[1];
rz(-2.3122687) q[1];
sx q[1];
rz(1.7620128) q[1];
x q[2];
rz(0.44646397) q[3];
sx q[3];
rz(-2.3013448) q[3];
sx q[3];
rz(-1.3069168) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(0.77439857) q[2];
sx q[2];
rz(-2.6937679) q[2];
sx q[2];
rz(1.1773342) q[2];
rz(2.2199953) q[3];
sx q[3];
rz(-2.2955743) q[3];
sx q[3];
rz(0.53068501) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.069799066) q[0];
sx q[0];
rz(-1.8853747) q[0];
sx q[0];
rz(1.0569093) q[0];
rz(-0.22843703) q[1];
sx q[1];
rz(-2.3616796) q[1];
sx q[1];
rz(2.8905919) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.042648249) q[0];
sx q[0];
rz(-1.3146719) q[0];
sx q[0];
rz(-2.8471208) q[0];
rz(0.86871882) q[2];
sx q[2];
rz(-1.0660604) q[2];
sx q[2];
rz(0.87930337) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.9343379) q[1];
sx q[1];
rz(-1.435507) q[1];
sx q[1];
rz(-2.9771718) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.9422705) q[3];
sx q[3];
rz(-2.0544996) q[3];
sx q[3];
rz(2.7886645) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-1.2061657) q[2];
sx q[2];
rz(-2.3212104) q[2];
sx q[2];
rz(2.7247735) q[2];
rz(2.5721512) q[3];
sx q[3];
rz(-2.0408401) q[3];
sx q[3];
rz(0.69863629) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.57132974) q[0];
sx q[0];
rz(-1.2705734) q[0];
sx q[0];
rz(0.50462333) q[0];
rz(-0.2568256) q[1];
sx q[1];
rz(-0.80373126) q[1];
sx q[1];
rz(-1.1508734) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9872091) q[0];
sx q[0];
rz(-2.0618274) q[0];
sx q[0];
rz(0.23316292) q[0];
rz(-1.1033642) q[2];
sx q[2];
rz(-2.4515477) q[2];
sx q[2];
rz(2.0375117) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.8650035) q[1];
sx q[1];
rz(-2.114608) q[1];
sx q[1];
rz(-2.0786042) q[1];
rz(-pi) q[2];
rz(-2.7959149) q[3];
sx q[3];
rz(-0.44000235) q[3];
sx q[3];
rz(-1.0087165) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.63200942) q[2];
sx q[2];
rz(-2.7183967) q[2];
sx q[2];
rz(-0.43935856) q[2];
rz(-1.1257233) q[3];
sx q[3];
rz(-1.6042387) q[3];
sx q[3];
rz(1.2996947) q[3];
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
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.87089649) q[0];
sx q[0];
rz(-0.82042158) q[0];
sx q[0];
rz(-2.6743555) q[0];
rz(-1.0386508) q[1];
sx q[1];
rz(-0.50913441) q[1];
sx q[1];
rz(-1.4899303) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2307546) q[0];
sx q[0];
rz(-2.0607134) q[0];
sx q[0];
rz(-1.3719489) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.9949399) q[2];
sx q[2];
rz(-1.6645558) q[2];
sx q[2];
rz(-2.3921426) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(2.6789602) q[1];
sx q[1];
rz(-2.3255146) q[1];
sx q[1];
rz(3.1051293) q[1];
x q[2];
rz(-0.64405264) q[3];
sx q[3];
rz(-2.4706512) q[3];
sx q[3];
rz(1.3356127) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.44522875) q[2];
sx q[2];
rz(-2.7859521) q[2];
sx q[2];
rz(-1.5871555) q[2];
rz(0.98443952) q[3];
sx q[3];
rz(-2.2319904) q[3];
sx q[3];
rz(-1.7279153) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5928891) q[0];
sx q[0];
rz(-2.2991572) q[0];
sx q[0];
rz(-0.62826759) q[0];
rz(-0.20333044) q[1];
sx q[1];
rz(-1.1140099) q[1];
sx q[1];
rz(-2.5471953) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3056476) q[0];
sx q[0];
rz(-1.3653269) q[0];
sx q[0];
rz(2.1688658) q[0];
x q[1];
rz(-0.35697414) q[2];
sx q[2];
rz(-2.5737615) q[2];
sx q[2];
rz(-0.50895509) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-0.59692057) q[1];
sx q[1];
rz(-1.871875) q[1];
sx q[1];
rz(1.8543509) q[1];
rz(-pi) q[2];
x q[2];
rz(1.6420308) q[3];
sx q[3];
rz(-2.1826577) q[3];
sx q[3];
rz(-0.16779403) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(1.2079042) q[2];
sx q[2];
rz(-0.7242569) q[2];
sx q[2];
rz(-0.17624632) q[2];
rz(-1.3147973) q[3];
sx q[3];
rz(-2.3254471) q[3];
sx q[3];
rz(0.76752457) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9888196) q[0];
sx q[0];
rz(-1.4403227) q[0];
sx q[0];
rz(1.3829917) q[0];
rz(3.026961) q[1];
sx q[1];
rz(-0.75420598) q[1];
sx q[1];
rz(2.4892714) q[1];
rz(-0.23723142) q[2];
sx q[2];
rz(-1.6085515) q[2];
sx q[2];
rz(-1.5534437) q[2];
rz(0.92656814) q[3];
sx q[3];
rz(-0.59048548) q[3];
sx q[3];
rz(0.98248972) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
