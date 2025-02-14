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
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7312429) q[0];
sx q[0];
rz(-2.2655537) q[0];
sx q[0];
rz(1.724051) q[0];
rz(-pi) q[1];
rz(1.4058787) q[2];
sx q[2];
rz(-1.5497297) q[2];
sx q[2];
rz(-1.5677468) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-0.47157745) q[1];
sx q[1];
rz(-1.0395998) q[1];
sx q[1];
rz(1.0047381) q[1];
rz(-pi) q[2];
rz(0.87907531) q[3];
sx q[3];
rz(-1.4180776) q[3];
sx q[3];
rz(-0.16593753) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(2.2012653) q[2];
sx q[2];
rz(-1.6028812) q[2];
sx q[2];
rz(0.24844696) q[2];
rz(1.1343608) q[3];
sx q[3];
rz(-3.0076707) q[3];
sx q[3];
rz(1.1321446) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.092875384) q[0];
sx q[0];
rz(-2.5889914) q[0];
sx q[0];
rz(-1.5930814) q[0];
rz(0.50367194) q[1];
sx q[1];
rz(-1.9182938) q[1];
sx q[1];
rz(0.37685397) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5232613) q[0];
sx q[0];
rz(-2.0428162) q[0];
sx q[0];
rz(1.5923772) q[0];
x q[1];
rz(-2.4433369) q[2];
sx q[2];
rz(-0.94329903) q[2];
sx q[2];
rz(2.1098532) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(0.26793626) q[1];
sx q[1];
rz(-2.5940478) q[1];
sx q[1];
rz(2.646562) q[1];
rz(-pi) q[2];
rz(-1.6474071) q[3];
sx q[3];
rz(-1.239068) q[3];
sx q[3];
rz(2.2463617) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-2.5571931) q[2];
sx q[2];
rz(-0.69584766) q[2];
sx q[2];
rz(2.2759571) q[2];
rz(1.4731167) q[3];
sx q[3];
rz(-0.98554635) q[3];
sx q[3];
rz(-0.98751155) q[3];
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
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5675548) q[0];
sx q[0];
rz(-1.2873298) q[0];
sx q[0];
rz(-1.0512742) q[0];
rz(-1.4569262) q[1];
sx q[1];
rz(-1.3070062) q[1];
sx q[1];
rz(-1.5164703) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2896958) q[0];
sx q[0];
rz(-2.3075884) q[0];
sx q[0];
rz(-1.9547321) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.3651354) q[2];
sx q[2];
rz(-0.8469905) q[2];
sx q[2];
rz(3.0645097) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(1.9439924) q[1];
sx q[1];
rz(-1.9514582) q[1];
sx q[1];
rz(1.0684816) q[1];
rz(-pi) q[2];
rz(-2.4184035) q[3];
sx q[3];
rz(-3.0602877) q[3];
sx q[3];
rz(-1.7618084) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-0.39525825) q[2];
sx q[2];
rz(-1.0552152) q[2];
sx q[2];
rz(0.22221097) q[2];
rz(-0.5018417) q[3];
sx q[3];
rz(-1.7343438) q[3];
sx q[3];
rz(-0.80880222) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.89268452) q[0];
sx q[0];
rz(-0.48766708) q[0];
sx q[0];
rz(-1.1647613) q[0];
rz(0.89272967) q[1];
sx q[1];
rz(-1.3803394) q[1];
sx q[1];
rz(-2.8597615) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.079275) q[0];
sx q[0];
rz(-1.6910416) q[0];
sx q[0];
rz(-2.6577302) q[0];
rz(-pi) q[1];
rz(-0.71433432) q[2];
sx q[2];
rz(-1.1118982) q[2];
sx q[2];
rz(0.66633505) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.6960245) q[1];
sx q[1];
rz(-0.17754517) q[1];
sx q[1];
rz(1.976718) q[1];
rz(-1.124566) q[3];
sx q[3];
rz(-1.5939624) q[3];
sx q[3];
rz(1.5678101) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-1.07936) q[2];
sx q[2];
rz(-0.56228176) q[2];
sx q[2];
rz(0.53323659) q[2];
rz(1.0821292) q[3];
sx q[3];
rz(-0.60492587) q[3];
sx q[3];
rz(2.3281039) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.37394062) q[0];
sx q[0];
rz(-2.7426608) q[0];
sx q[0];
rz(1.8178513) q[0];
rz(2.3228877) q[1];
sx q[1];
rz(-2.3981514) q[1];
sx q[1];
rz(-1.0680107) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6477953) q[0];
sx q[0];
rz(-1.6500705) q[0];
sx q[0];
rz(0.1597516) q[0];
rz(-pi) q[1];
rz(-0.83663656) q[2];
sx q[2];
rz(-1.0435487) q[2];
sx q[2];
rz(-0.0076779445) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.45324486) q[1];
sx q[1];
rz(-1.6736341) q[1];
sx q[1];
rz(-2.4842841) q[1];
rz(2.976494) q[3];
sx q[3];
rz(-2.0378651) q[3];
sx q[3];
rz(1.1822753) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-1.1800804) q[2];
sx q[2];
rz(-1.1137806) q[2];
sx q[2];
rz(2.2799344) q[2];
rz(1.0263475) q[3];
sx q[3];
rz(-1.15871) q[3];
sx q[3];
rz(2.0238743) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3468129) q[0];
sx q[0];
rz(-2.3220799) q[0];
sx q[0];
rz(-0.89282194) q[0];
rz(-2.9606441) q[1];
sx q[1];
rz(-0.67829689) q[1];
sx q[1];
rz(-0.13042626) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3149726) q[0];
sx q[0];
rz(-2.9493647) q[0];
sx q[0];
rz(0.97127025) q[0];
x q[1];
rz(2.4685235) q[2];
sx q[2];
rz(-2.4166738) q[2];
sx q[2];
rz(-0.47436213) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.0196258) q[1];
sx q[1];
rz(-1.6000012) q[1];
sx q[1];
rz(1.8080416) q[1];
rz(-0.78108861) q[3];
sx q[3];
rz(-1.4295331) q[3];
sx q[3];
rz(-3.0716592) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.22455198) q[2];
sx q[2];
rz(-1.2265393) q[2];
sx q[2];
rz(-0.87635931) q[2];
rz(1.2419491) q[3];
sx q[3];
rz(-0.94049898) q[3];
sx q[3];
rz(1.0103286) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4001615) q[0];
sx q[0];
rz(-3.042996) q[0];
sx q[0];
rz(-0.36369351) q[0];
rz(-0.74186507) q[1];
sx q[1];
rz(-2.1153085) q[1];
sx q[1];
rz(0.40282869) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8217709) q[0];
sx q[0];
rz(-1.9812225) q[0];
sx q[0];
rz(-2.8794146) q[0];
x q[1];
rz(2.7637787) q[2];
sx q[2];
rz(-1.1973518) q[2];
sx q[2];
rz(0.4682954) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(0.070870474) q[1];
sx q[1];
rz(-1.6477343) q[1];
sx q[1];
rz(-1.6926195) q[1];
rz(1.2996583) q[3];
sx q[3];
rz(-0.95356546) q[3];
sx q[3];
rz(-0.24880508) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(1.7824629) q[2];
sx q[2];
rz(-2.7733347) q[2];
sx q[2];
rz(-1.5943607) q[2];
rz(-2.3063229) q[3];
sx q[3];
rz(-0.95649496) q[3];
sx q[3];
rz(-0.48867759) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7816724) q[0];
sx q[0];
rz(-1.3363573) q[0];
sx q[0];
rz(-0.51505995) q[0];
rz(1.4393648) q[1];
sx q[1];
rz(-1.2934338) q[1];
sx q[1];
rz(-2.7424367) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.24267749) q[0];
sx q[0];
rz(-1.5726046) q[0];
sx q[0];
rz(3.010672) q[0];
x q[1];
rz(-0.24243124) q[2];
sx q[2];
rz(-2.2848086) q[2];
sx q[2];
rz(-2.8664608) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.29447435) q[1];
sx q[1];
rz(-0.26853466) q[1];
sx q[1];
rz(1.8163791) q[1];
rz(-pi) q[2];
rz(0.065654556) q[3];
sx q[3];
rz(-1.5720417) q[3];
sx q[3];
rz(-0.37855442) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-1.9814375) q[2];
sx q[2];
rz(-1.258054) q[2];
sx q[2];
rz(1.0464279) q[2];
rz(0.6662755) q[3];
sx q[3];
rz(-0.25686887) q[3];
sx q[3];
rz(-1.0506786) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
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
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.94205725) q[0];
sx q[0];
rz(-0.10534795) q[0];
sx q[0];
rz(0.60607213) q[0];
rz(-2.3145158) q[1];
sx q[1];
rz(-1.3304293) q[1];
sx q[1];
rz(-1.3333295) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5510941) q[0];
sx q[0];
rz(-2.7685389) q[0];
sx q[0];
rz(-2.8540552) q[0];
x q[1];
rz(-0.74335783) q[2];
sx q[2];
rz(-1.2326692) q[2];
sx q[2];
rz(1.0184792) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(3.044802) q[1];
sx q[1];
rz(-1.6875182) q[1];
sx q[1];
rz(2.2404284) q[1];
rz(-pi) q[2];
rz(1.656084) q[3];
sx q[3];
rz(-1.1138289) q[3];
sx q[3];
rz(-2.5449139) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(1.0102319) q[2];
sx q[2];
rz(-2.2183245) q[2];
sx q[2];
rz(2.7566747) q[2];
rz(-0.63079232) q[3];
sx q[3];
rz(-1.2814949) q[3];
sx q[3];
rz(1.3425286) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
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
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.25403062) q[0];
sx q[0];
rz(-1.3986724) q[0];
sx q[0];
rz(-1.4400462) q[0];
rz(0.74132672) q[1];
sx q[1];
rz(-1.4011551) q[1];
sx q[1];
rz(0.25873605) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.16782204) q[0];
sx q[0];
rz(-2.8956465) q[0];
sx q[0];
rz(1.5542547) q[0];
x q[1];
rz(-1.0526686) q[2];
sx q[2];
rz(-1.8473704) q[2];
sx q[2];
rz(2.5605048) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(2.1867793) q[1];
sx q[1];
rz(-2.3329321) q[1];
sx q[1];
rz(2.0524128) q[1];
rz(-pi) q[2];
rz(0.40269884) q[3];
sx q[3];
rz(-0.81953632) q[3];
sx q[3];
rz(2.0961026) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-2.4197293) q[2];
sx q[2];
rz(-1.1121007) q[2];
sx q[2];
rz(1.8224243) q[2];
rz(-2.3223274) q[3];
sx q[3];
rz(-2.0115439) q[3];
sx q[3];
rz(-0.54660249) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4010314) q[0];
sx q[0];
rz(-0.87420976) q[0];
sx q[0];
rz(2.5921205) q[0];
rz(1.4389379) q[1];
sx q[1];
rz(-0.66822744) q[1];
sx q[1];
rz(-2.6034036) q[1];
rz(1.5835252) q[2];
sx q[2];
rz(-0.71011484) q[2];
sx q[2];
rz(-1.9488123) q[2];
rz(-0.2470369) q[3];
sx q[3];
rz(-1.023019) q[3];
sx q[3];
rz(-1.116397) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
