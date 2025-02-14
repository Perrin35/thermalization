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
rz(2.574805) q[0];
rz(1.2504638) q[1];
sx q[1];
rz(-1.992978) q[1];
sx q[1];
rz(2.221938) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.56432589) q[0];
sx q[0];
rz(-0.79123215) q[0];
sx q[0];
rz(1.3444232) q[0];
x q[1];
rz(-0.55962015) q[2];
sx q[2];
rz(-2.3336448) q[2];
sx q[2];
rz(1.1839649) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(0.39988841) q[1];
sx q[1];
rz(-0.9298069) q[1];
sx q[1];
rz(-1.5128894) q[1];
rz(-pi) q[2];
rz(1.4266495) q[3];
sx q[3];
rz(-1.0394154) q[3];
sx q[3];
rz(1.2856158) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(2.2016242) q[2];
sx q[2];
rz(-2.2087966) q[2];
sx q[2];
rz(2.0882108) q[2];
rz(2.4273704) q[3];
sx q[3];
rz(-0.37319365) q[3];
sx q[3];
rz(-1.2936973) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
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
rz(1.6222222) q[0];
sx q[0];
rz(-0.70835963) q[0];
sx q[0];
rz(0.57193065) q[0];
rz(1.2988623) q[1];
sx q[1];
rz(-2.3426901) q[1];
sx q[1];
rz(-0.78972185) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2970726) q[0];
sx q[0];
rz(-1.6668586) q[0];
sx q[0];
rz(1.8734832) q[0];
x q[1];
rz(2.4873494) q[2];
sx q[2];
rz(-1.4432608) q[2];
sx q[2];
rz(-2.5728284) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.8487723) q[1];
sx q[1];
rz(-2.1629984) q[1];
sx q[1];
rz(1.2315537) q[1];
rz(-0.67482194) q[3];
sx q[3];
rz(-1.1742419) q[3];
sx q[3];
rz(-3.0930072) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-2.2907437) q[2];
sx q[2];
rz(-2.3727543) q[2];
sx q[2];
rz(1.7791465) q[2];
rz(-0.29081523) q[3];
sx q[3];
rz(-1.3664061) q[3];
sx q[3];
rz(-3.0349558) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
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
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.13880759) q[0];
sx q[0];
rz(-1.0951575) q[0];
sx q[0];
rz(0.7005257) q[0];
rz(-0.48577148) q[1];
sx q[1];
rz(-2.3120717) q[1];
sx q[1];
rz(-0.22451678) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1248847) q[0];
sx q[0];
rz(-1.2205692) q[0];
sx q[0];
rz(1.4136397) q[0];
x q[1];
rz(3.0430138) q[2];
sx q[2];
rz(-0.99462992) q[2];
sx q[2];
rz(-0.88159943) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-0.77188801) q[1];
sx q[1];
rz(-1.4955031) q[1];
sx q[1];
rz(2.1545707) q[1];
rz(-pi) q[2];
x q[2];
rz(0.84214476) q[3];
sx q[3];
rz(-1.6557459) q[3];
sx q[3];
rz(-0.85499518) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.38480467) q[2];
sx q[2];
rz(-2.2497358) q[2];
sx q[2];
rz(3.0925114) q[2];
rz(-0.067042025) q[3];
sx q[3];
rz(-1.9837572) q[3];
sx q[3];
rz(2.4750347) q[3];
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
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9662358) q[0];
sx q[0];
rz(-0.55713621) q[0];
sx q[0];
rz(3.1224342) q[0];
rz(-1.3592023) q[1];
sx q[1];
rz(-1.7117585) q[1];
sx q[1];
rz(-0.0044936831) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.053872975) q[0];
sx q[0];
rz(-1.036265) q[0];
sx q[0];
rz(0.9854395) q[0];
rz(0.7257113) q[2];
sx q[2];
rz(-2.5134676) q[2];
sx q[2];
rz(-0.84838644) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(1.5792721) q[1];
sx q[1];
rz(-1.4021789) q[1];
sx q[1];
rz(-1.2071441) q[1];
x q[2];
rz(-1.6636623) q[3];
sx q[3];
rz(-1.3680653) q[3];
sx q[3];
rz(0.40846881) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-2.6898474) q[2];
sx q[2];
rz(-1.0901901) q[2];
sx q[2];
rz(-1.4534265) q[2];
rz(1.2545741) q[3];
sx q[3];
rz(-1.6202319) q[3];
sx q[3];
rz(-0.5622676) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3314064) q[0];
sx q[0];
rz(-1.2368546) q[0];
sx q[0];
rz(-1.1619262) q[0];
rz(-0.76238531) q[1];
sx q[1];
rz(-0.98680174) q[1];
sx q[1];
rz(0.42957482) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9305161) q[0];
sx q[0];
rz(-1.9787702) q[0];
sx q[0];
rz(3.1019449) q[0];
rz(-pi) q[1];
rz(2.863071) q[2];
sx q[2];
rz(-2.3633011) q[2];
sx q[2];
rz(0.24519224) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.9813198) q[1];
sx q[1];
rz(-1.2339051) q[1];
sx q[1];
rz(-0.87967061) q[1];
rz(2.0031092) q[3];
sx q[3];
rz(-0.70037445) q[3];
sx q[3];
rz(-2.5770503) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-0.21542159) q[2];
sx q[2];
rz(-1.8502219) q[2];
sx q[2];
rz(1.311709) q[2];
rz(-0.46071509) q[3];
sx q[3];
rz(-2.1381133) q[3];
sx q[3];
rz(-0.13319143) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8120414) q[0];
sx q[0];
rz(-2.9819745) q[0];
sx q[0];
rz(2.0352236) q[0];
rz(0.1144935) q[1];
sx q[1];
rz(-2.0888927) q[1];
sx q[1];
rz(0.053650275) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6686443) q[0];
sx q[0];
rz(-1.8264174) q[0];
sx q[0];
rz(0.0059519569) q[0];
rz(0.80751632) q[2];
sx q[2];
rz(-1.6431333) q[2];
sx q[2];
rz(1.4405719) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-2.304639) q[1];
sx q[1];
rz(-2.0817679) q[1];
sx q[1];
rz(-3.0797019) q[1];
x q[2];
rz(-2.3785081) q[3];
sx q[3];
rz(-0.79417919) q[3];
sx q[3];
rz(-1.1556243) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(1.2758808) q[2];
sx q[2];
rz(-1.2578473) q[2];
sx q[2];
rz(-2.6344521) q[2];
rz(-1.3745314) q[3];
sx q[3];
rz(-1.6009228) q[3];
sx q[3];
rz(1.1486294) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.62992612) q[0];
sx q[0];
rz(-2.0138795) q[0];
sx q[0];
rz(0.041286904) q[0];
rz(2.4033026) q[1];
sx q[1];
rz(-1.2246882) q[1];
sx q[1];
rz(-0.98947492) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6302297) q[0];
sx q[0];
rz(-2.0718899) q[0];
sx q[0];
rz(2.5998235) q[0];
rz(-0.35571738) q[2];
sx q[2];
rz(-1.748744) q[2];
sx q[2];
rz(-1.8989814) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.38997) q[1];
sx q[1];
rz(-0.97725673) q[1];
sx q[1];
rz(-2.4291527) q[1];
rz(-pi) q[2];
x q[2];
rz(2.2253898) q[3];
sx q[3];
rz(-1.8713017) q[3];
sx q[3];
rz(-0.88778561) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-1.5422633) q[2];
sx q[2];
rz(-1.0487391) q[2];
sx q[2];
rz(0.85912022) q[2];
rz(3.1298992) q[3];
sx q[3];
rz(-1.6418567) q[3];
sx q[3];
rz(-3.0288127) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3149253) q[0];
sx q[0];
rz(-2.5496917) q[0];
sx q[0];
rz(-2.9097606) q[0];
rz(2.5595317) q[1];
sx q[1];
rz(-1.8128017) q[1];
sx q[1];
rz(1.9225165) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9052637) q[0];
sx q[0];
rz(-1.4538611) q[0];
sx q[0];
rz(0.56044436) q[0];
rz(-pi) q[1];
x q[1];
rz(0.14149551) q[2];
sx q[2];
rz(-2.7853051) q[2];
sx q[2];
rz(0.74559272) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.2197709) q[1];
sx q[1];
rz(-1.3857462) q[1];
sx q[1];
rz(-0.15404032) q[1];
rz(-2.3858504) q[3];
sx q[3];
rz(-1.1583987) q[3];
sx q[3];
rz(-0.88471593) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.9570738) q[2];
sx q[2];
rz(-1.3513214) q[2];
sx q[2];
rz(-0.78682023) q[2];
rz(1.6400853) q[3];
sx q[3];
rz(-2.7101176) q[3];
sx q[3];
rz(-1.4260346) q[3];
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
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.96104923) q[0];
sx q[0];
rz(-1.3732055) q[0];
sx q[0];
rz(-0.94026646) q[0];
rz(-0.44479784) q[1];
sx q[1];
rz(-0.85591379) q[1];
sx q[1];
rz(0.14642265) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.0068374182) q[0];
sx q[0];
rz(-1.869258) q[0];
sx q[0];
rz(1.439007) q[0];
rz(-pi) q[1];
x q[1];
rz(0.27250473) q[2];
sx q[2];
rz(-1.3256058) q[2];
sx q[2];
rz(1.6134855) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-0.89288974) q[1];
sx q[1];
rz(-1.8161402) q[1];
sx q[1];
rz(2.4470083) q[1];
rz(-pi) q[2];
rz(3.0412263) q[3];
sx q[3];
rz(-2.3477738) q[3];
sx q[3];
rz(2.9507) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-3.0545097) q[2];
sx q[2];
rz(-1.4887709) q[2];
sx q[2];
rz(-2.3040237) q[2];
rz(-0.93291035) q[3];
sx q[3];
rz(-0.98740238) q[3];
sx q[3];
rz(-2.3605409) q[3];
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
x q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8785716) q[0];
sx q[0];
rz(-1.2282547) q[0];
sx q[0];
rz(1.7735057) q[0];
rz(-3.0655762) q[1];
sx q[1];
rz(-0.58034211) q[1];
sx q[1];
rz(1.1200294) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.0033542643) q[0];
sx q[0];
rz(-1.3068259) q[0];
sx q[0];
rz(1.2498943) q[0];
rz(0.70856382) q[2];
sx q[2];
rz(-1.8052793) q[2];
sx q[2];
rz(1.3402779) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(1.315584) q[1];
sx q[1];
rz(-0.30696973) q[1];
sx q[1];
rz(1.2493285) q[1];
rz(-pi) q[2];
rz(-1.7342274) q[3];
sx q[3];
rz(-1.4686958) q[3];
sx q[3];
rz(-0.40615505) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.7398305) q[2];
sx q[2];
rz(-1.4515452) q[2];
sx q[2];
rz(2.8945727) q[2];
rz(0.57010993) q[3];
sx q[3];
rz(-0.48738042) q[3];
sx q[3];
rz(-2.0232239) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0433255) q[0];
sx q[0];
rz(-1.387431) q[0];
sx q[0];
rz(-0.36547216) q[0];
rz(-0.098516057) q[1];
sx q[1];
rz(-2.5010074) q[1];
sx q[1];
rz(0.88190257) q[1];
rz(-3.1280853) q[2];
sx q[2];
rz(-1.1579222) q[2];
sx q[2];
rz(2.9029878) q[2];
rz(-2.5539342) q[3];
sx q[3];
rz(-2.6937204) q[3];
sx q[3];
rz(3.1093521) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
