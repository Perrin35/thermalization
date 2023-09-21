OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-1.9010889) q[0];
sx q[0];
rz(-0.17740372) q[0];
sx q[0];
rz(-2.0071964) q[0];
rz(-1.9534684) q[1];
sx q[1];
rz(-1.0367353) q[1];
sx q[1];
rz(-2.477975) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.21270574) q[0];
sx q[0];
rz(-2.0318673) q[0];
sx q[0];
rz(1.5972208) q[0];
x q[1];
rz(-1.9248336) q[2];
sx q[2];
rz(-0.82497043) q[2];
sx q[2];
rz(2.0113457) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-2.8780958) q[1];
sx q[1];
rz(-0.43259183) q[1];
sx q[1];
rz(2.2357975) q[1];
rz(0.10856467) q[3];
sx q[3];
rz(-1.3285471) q[3];
sx q[3];
rz(0.066075174) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.87876451) q[2];
sx q[2];
rz(-0.44439134) q[2];
sx q[2];
rz(3.0905241) q[2];
rz(-2.5845394) q[3];
sx q[3];
rz(-0.80018187) q[3];
sx q[3];
rz(-1.5867656) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.58650815) q[0];
sx q[0];
rz(-2.3095135) q[0];
sx q[0];
rz(-2.5449975) q[0];
rz(2.3157628) q[1];
sx q[1];
rz(-1.4412216) q[1];
sx q[1];
rz(-1.2260431) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.37218371) q[0];
sx q[0];
rz(-0.49159494) q[0];
sx q[0];
rz(-2.7093637) q[0];
rz(-pi) q[1];
rz(2.3833582) q[2];
sx q[2];
rz(-2.1633534) q[2];
sx q[2];
rz(-0.50298467) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-3.1184428) q[1];
sx q[1];
rz(-0.21986248) q[1];
sx q[1];
rz(-1.5688194) q[1];
rz(-pi) q[2];
x q[2];
rz(-3.1321208) q[3];
sx q[3];
rz(-2.1433899) q[3];
sx q[3];
rz(-0.026281683) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(0.79919672) q[2];
sx q[2];
rz(-1.1725972) q[2];
sx q[2];
rz(-2.8105695) q[2];
rz(0.80667574) q[3];
sx q[3];
rz(-1.9162063) q[3];
sx q[3];
rz(-1.4276918) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
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
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1598635) q[0];
sx q[0];
rz(-1.9112497) q[0];
sx q[0];
rz(1.8925517) q[0];
rz(-3.0535835) q[1];
sx q[1];
rz(-2.0188589) q[1];
sx q[1];
rz(1.0294611) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2217076) q[0];
sx q[0];
rz(-1.3516597) q[0];
sx q[0];
rz(-0.94421454) q[0];
rz(-pi) q[1];
x q[1];
rz(2.4918409) q[2];
sx q[2];
rz(-1.8012815) q[2];
sx q[2];
rz(2.8284555) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(0.50476915) q[1];
sx q[1];
rz(-1.4154589) q[1];
sx q[1];
rz(1.4298646) q[1];
rz(-pi) q[2];
rz(0.73996468) q[3];
sx q[3];
rz(-1.4810586) q[3];
sx q[3];
rz(1.9268074) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-3.133698) q[2];
sx q[2];
rz(-1.4174613) q[2];
sx q[2];
rz(2.6339445) q[2];
rz(-1.3890022) q[3];
sx q[3];
rz(-0.32998431) q[3];
sx q[3];
rz(1.1631789) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4841109) q[0];
sx q[0];
rz(-2.8634475) q[0];
sx q[0];
rz(1.5959651) q[0];
rz(1.0428628) q[1];
sx q[1];
rz(-1.9680126) q[1];
sx q[1];
rz(1.5159336) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3665873) q[0];
sx q[0];
rz(-2.2457652) q[0];
sx q[0];
rz(-1.3875899) q[0];
rz(2.1298725) q[2];
sx q[2];
rz(-1.498073) q[2];
sx q[2];
rz(2.5474472) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-3.1277395) q[1];
sx q[1];
rz(-1.2989559) q[1];
sx q[1];
rz(-0.66687648) q[1];
rz(-pi) q[2];
x q[2];
rz(1.7901778) q[3];
sx q[3];
rz(-0.99321584) q[3];
sx q[3];
rz(2.1882309) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.76379124) q[2];
sx q[2];
rz(-0.8447454) q[2];
sx q[2];
rz(-1.2949004) q[2];
rz(0.14136782) q[3];
sx q[3];
rz(-0.54261345) q[3];
sx q[3];
rz(-0.98658371) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7871053) q[0];
sx q[0];
rz(-1.961345) q[0];
sx q[0];
rz(1.6217344) q[0];
rz(-0.63201085) q[1];
sx q[1];
rz(-2.4193587) q[1];
sx q[1];
rz(0.89486665) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.97840727) q[0];
sx q[0];
rz(-2.1911231) q[0];
sx q[0];
rz(2.083367) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.52389) q[2];
sx q[2];
rz(-2.3139944) q[2];
sx q[2];
rz(-2.655381) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(1.4836854) q[1];
sx q[1];
rz(-2.2135418) q[1];
sx q[1];
rz(0.55773463) q[1];
x q[2];
rz(-0.28646333) q[3];
sx q[3];
rz(-1.5634007) q[3];
sx q[3];
rz(1.8252107) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(0.52508369) q[2];
sx q[2];
rz(-0.36642763) q[2];
sx q[2];
rz(-1.1425225) q[2];
rz(2.3948495) q[3];
sx q[3];
rz(-1.2544422) q[3];
sx q[3];
rz(2.0660627) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.557945) q[0];
sx q[0];
rz(-0.32145158) q[0];
sx q[0];
rz(-1.3775795) q[0];
rz(-2.6691943) q[1];
sx q[1];
rz(-0.51858416) q[1];
sx q[1];
rz(-0.46498743) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8573479) q[0];
sx q[0];
rz(-1.9042958) q[0];
sx q[0];
rz(2.0302982) q[0];
rz(2.8000185) q[2];
sx q[2];
rz(-1.6401059) q[2];
sx q[2];
rz(2.9937033) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-1.2496693) q[1];
sx q[1];
rz(-1.6185456) q[1];
sx q[1];
rz(0.93530099) q[1];
rz(0.70907866) q[3];
sx q[3];
rz(-1.7311829) q[3];
sx q[3];
rz(-1.3080314) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-0.49089367) q[2];
sx q[2];
rz(-1.2630254) q[2];
sx q[2];
rz(-0.4450376) q[2];
rz(-0.93368357) q[3];
sx q[3];
rz(-1.4327587) q[3];
sx q[3];
rz(0.26708189) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.12748195) q[0];
sx q[0];
rz(-1.5662136) q[0];
sx q[0];
rz(-1.7215464) q[0];
rz(0.02380112) q[1];
sx q[1];
rz(-0.61444608) q[1];
sx q[1];
rz(-0.15596095) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.10931817) q[0];
sx q[0];
rz(-0.30800691) q[0];
sx q[0];
rz(0.59662915) q[0];
rz(0.51257001) q[2];
sx q[2];
rz(-1.9480431) q[2];
sx q[2];
rz(0.23114983) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-1.5739463) q[1];
sx q[1];
rz(-1.624755) q[1];
sx q[1];
rz(1.3957363) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.53653036) q[3];
sx q[3];
rz(-1.1324258) q[3];
sx q[3];
rz(-2.8560672) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(0.069313958) q[2];
sx q[2];
rz(-2.5645655) q[2];
sx q[2];
rz(-1.0726661) q[2];
rz(0.32564751) q[3];
sx q[3];
rz(-1.1762534) q[3];
sx q[3];
rz(1.5163039) q[3];
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
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1919365) q[0];
sx q[0];
rz(-0.40238109) q[0];
sx q[0];
rz(0.33777133) q[0];
rz(-1.0900963) q[1];
sx q[1];
rz(-0.97507674) q[1];
sx q[1];
rz(2.8930194) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0524806) q[0];
sx q[0];
rz(-2.5626474) q[0];
sx q[0];
rz(2.6738033) q[0];
rz(-2.3209004) q[2];
sx q[2];
rz(-2.4889915) q[2];
sx q[2];
rz(-1.5479659) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(3.0174745) q[1];
sx q[1];
rz(-1.0867456) q[1];
sx q[1];
rz(-2.2494621) q[1];
rz(-pi) q[2];
rz(-1.7406086) q[3];
sx q[3];
rz(-1.4949993) q[3];
sx q[3];
rz(-2.4236987) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(1.8481855) q[2];
sx q[2];
rz(-2.6082787) q[2];
sx q[2];
rz(-0.08671134) q[2];
rz(-2.6596206) q[3];
sx q[3];
rz(-0.50589365) q[3];
sx q[3];
rz(1.988525) q[3];
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
rz(pi/2) q[0];
rz(-pi) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6329353) q[0];
sx q[0];
rz(-1.0077227) q[0];
sx q[0];
rz(2.8588262) q[0];
rz(-0.70156082) q[1];
sx q[1];
rz(-0.82156721) q[1];
sx q[1];
rz(1.3185906) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0193664) q[0];
sx q[0];
rz(-1.1302233) q[0];
sx q[0];
rz(-0.088935436) q[0];
rz(-pi) q[1];
rz(-0.96005) q[2];
sx q[2];
rz(-2.1337482) q[2];
sx q[2];
rz(0.27573953) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(0.56837592) q[1];
sx q[1];
rz(-1.425256) q[1];
sx q[1];
rz(0.88370609) q[1];
x q[2];
rz(-1.1845469) q[3];
sx q[3];
rz(-0.85856122) q[3];
sx q[3];
rz(-1.4854747) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(2.7302154) q[2];
sx q[2];
rz(-2.8432379) q[2];
sx q[2];
rz(2.7424157) q[2];
rz(-0.88360751) q[3];
sx q[3];
rz(-1.6688321) q[3];
sx q[3];
rz(1.2214899) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3783962) q[0];
sx q[0];
rz(-2.3662687) q[0];
sx q[0];
rz(-1.5989074) q[0];
rz(-2.0762766) q[1];
sx q[1];
rz(-0.97384802) q[1];
sx q[1];
rz(1.261196) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6529918) q[0];
sx q[0];
rz(-1.4654667) q[0];
sx q[0];
rz(-0.41098849) q[0];
rz(1.5930428) q[2];
sx q[2];
rz(-0.97264475) q[2];
sx q[2];
rz(-1.0215789) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(2.2640472) q[1];
sx q[1];
rz(-2.0685158) q[1];
sx q[1];
rz(1.4977786) q[1];
rz(-2.4576549) q[3];
sx q[3];
rz(-1.4468907) q[3];
sx q[3];
rz(-2.370196) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(2.5836872) q[2];
sx q[2];
rz(-2.5186899) q[2];
sx q[2];
rz(2.5718001) q[2];
rz(-1.9231046) q[3];
sx q[3];
rz(-2.4945538) q[3];
sx q[3];
rz(-0.56263721) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
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
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8284843) q[0];
sx q[0];
rz(-2.2611571) q[0];
sx q[0];
rz(-1.8631998) q[0];
rz(2.5333511) q[1];
sx q[1];
rz(-0.47641644) q[1];
sx q[1];
rz(0.48412916) q[1];
rz(-3.1153395) q[2];
sx q[2];
rz(-2.142341) q[2];
sx q[2];
rz(-0.061948902) q[2];
rz(-2.5504997) q[3];
sx q[3];
rz(-1.686284) q[3];
sx q[3];
rz(-1.4808663) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];