OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-0.91575032) q[0];
sx q[0];
rz(3.1728035) q[0];
sx q[0];
rz(6.7682545) q[0];
rz(-2.354061) q[1];
sx q[1];
rz(-2.1252316) q[1];
sx q[1];
rz(0.41419849) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.83672268) q[0];
sx q[0];
rz(-0.49855907) q[0];
sx q[0];
rz(-0.81122938) q[0];
rz(-1.0010927) q[2];
sx q[2];
rz(-1.4909407) q[2];
sx q[2];
rz(-0.18730883) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.5536082) q[1];
sx q[1];
rz(-1.6349287) q[1];
sx q[1];
rz(-1.1448121) q[1];
rz(-pi) q[2];
rz(2.7798153) q[3];
sx q[3];
rz(-2.8594115) q[3];
sx q[3];
rz(-0.85229814) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(1.9238613) q[2];
sx q[2];
rz(-1.8863181) q[2];
sx q[2];
rz(-0.031575354) q[2];
rz(-1.2565553) q[3];
sx q[3];
rz(-2.6387408) q[3];
sx q[3];
rz(2.6149635) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8388222) q[0];
sx q[0];
rz(-1.4571723) q[0];
sx q[0];
rz(-2.9717428) q[0];
rz(-0.70392144) q[1];
sx q[1];
rz(-2.070065) q[1];
sx q[1];
rz(0.53952113) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8115494) q[0];
sx q[0];
rz(-1.1216315) q[0];
sx q[0];
rz(3.0898068) q[0];
x q[1];
rz(1.9794481) q[2];
sx q[2];
rz(-2.2505629) q[2];
sx q[2];
rz(-2.0603927) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.62944618) q[1];
sx q[1];
rz(-2.0680032) q[1];
sx q[1];
rz(-2.9260103) q[1];
rz(2.2803454) q[3];
sx q[3];
rz(-1.1303139) q[3];
sx q[3];
rz(-1.7064106) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-2.4743621) q[2];
sx q[2];
rz(-1.903406) q[2];
sx q[2];
rz(-2.898522) q[2];
rz(2.4754751) q[3];
sx q[3];
rz(-0.56454286) q[3];
sx q[3];
rz(-1.8977785) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.77984017) q[0];
sx q[0];
rz(-3.0267974) q[0];
sx q[0];
rz(0.4483805) q[0];
rz(1.7547296) q[1];
sx q[1];
rz(-1.9877537) q[1];
sx q[1];
rz(2.8853436) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6305144) q[0];
sx q[0];
rz(-0.96195463) q[0];
sx q[0];
rz(1.2468673) q[0];
x q[1];
rz(2.0492378) q[2];
sx q[2];
rz(-1.1698327) q[2];
sx q[2];
rz(1.5329597) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-0.33830723) q[1];
sx q[1];
rz(-0.76122621) q[1];
sx q[1];
rz(1.085698) q[1];
rz(-pi) q[2];
rz(0.20817169) q[3];
sx q[3];
rz(-2.9609207) q[3];
sx q[3];
rz(-0.65073035) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(1.3391352) q[2];
sx q[2];
rz(-2.4352303) q[2];
sx q[2];
rz(2.5668872) q[2];
rz(1.7859219) q[3];
sx q[3];
rz(-1.971743) q[3];
sx q[3];
rz(1.2333966) q[3];
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
x q[0];
rz(pi/2) q[0];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9280076) q[0];
sx q[0];
rz(-1.428823) q[0];
sx q[0];
rz(-2.8821049) q[0];
rz(1.150594) q[1];
sx q[1];
rz(-1.78803) q[1];
sx q[1];
rz(2.4096699) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4703092) q[0];
sx q[0];
rz(-1.0571612) q[0];
sx q[0];
rz(-2.1156103) q[0];
x q[1];
rz(0.71728431) q[2];
sx q[2];
rz(-1.4826164) q[2];
sx q[2];
rz(0.35536534) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(0.54938984) q[1];
sx q[1];
rz(-1.376774) q[1];
sx q[1];
rz(-0.26280304) q[1];
rz(2.6287574) q[3];
sx q[3];
rz(-0.25179112) q[3];
sx q[3];
rz(2.4076715) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(1.4449473) q[2];
sx q[2];
rz(-1.2735294) q[2];
sx q[2];
rz(-0.15110061) q[2];
rz(0.54667306) q[3];
sx q[3];
rz(-2.0918545) q[3];
sx q[3];
rz(2.7643519) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.54995173) q[0];
sx q[0];
rz(-2.0786091) q[0];
sx q[0];
rz(-1.8792101) q[0];
rz(-1.4683912) q[1];
sx q[1];
rz(-0.60931283) q[1];
sx q[1];
rz(2.343822) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7933465) q[0];
sx q[0];
rz(-3.0214546) q[0];
sx q[0];
rz(-1.5500463) q[0];
rz(-2.8380978) q[2];
sx q[2];
rz(-1.7804838) q[2];
sx q[2];
rz(-1.4022624) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-1.5835727) q[1];
sx q[1];
rz(-0.94545525) q[1];
sx q[1];
rz(-3.0261092) q[1];
rz(-pi) q[2];
rz(2.6990715) q[3];
sx q[3];
rz(-1.8707152) q[3];
sx q[3];
rz(2.4181441) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(2.795934) q[2];
sx q[2];
rz(-2.5107333) q[2];
sx q[2];
rz(-2.8395555) q[2];
rz(-1.9942412) q[3];
sx q[3];
rz(-1.4619504) q[3];
sx q[3];
rz(0.54774493) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7323332) q[0];
sx q[0];
rz(-3.0497666) q[0];
sx q[0];
rz(-1.1557895) q[0];
rz(2.0571158) q[1];
sx q[1];
rz(-0.98025727) q[1];
sx q[1];
rz(0.070080431) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.72567155) q[0];
sx q[0];
rz(-1.7416746) q[0];
sx q[0];
rz(0.15539774) q[0];
rz(-pi) q[1];
rz(-2.5885133) q[2];
sx q[2];
rz(-0.86420176) q[2];
sx q[2];
rz(-1.7318219) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-0.49680432) q[1];
sx q[1];
rz(-1.4092688) q[1];
sx q[1];
rz(0.5715538) q[1];
rz(-2.0134986) q[3];
sx q[3];
rz(-3.0590995) q[3];
sx q[3];
rz(0.024162956) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-0.30248102) q[2];
sx q[2];
rz(-1.182686) q[2];
sx q[2];
rz(-2.5202259) q[2];
rz(-1.7012043) q[3];
sx q[3];
rz(-0.50783235) q[3];
sx q[3];
rz(-0.5733718) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.086534111) q[0];
sx q[0];
rz(-1.4165514) q[0];
sx q[0];
rz(0.85987464) q[0];
rz(1.9372008) q[1];
sx q[1];
rz(-0.87221968) q[1];
sx q[1];
rz(3.133657) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.313109) q[0];
sx q[0];
rz(-2.5223753) q[0];
sx q[0];
rz(2.6909268) q[0];
x q[1];
rz(-1.0646348) q[2];
sx q[2];
rz(-0.75876615) q[2];
sx q[2];
rz(0.90492349) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(2.4850033) q[1];
sx q[1];
rz(-0.90066972) q[1];
sx q[1];
rz(0.90799241) q[1];
x q[2];
rz(-0.90585917) q[3];
sx q[3];
rz(-2.6568036) q[3];
sx q[3];
rz(-1.0672027) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.4456711) q[2];
sx q[2];
rz(-1.9182703) q[2];
sx q[2];
rz(2.725214) q[2];
rz(-1.773206) q[3];
sx q[3];
rz(-1.2973283) q[3];
sx q[3];
rz(2.2369475) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
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
rz(0.96173441) q[0];
sx q[0];
rz(-0.27357736) q[0];
sx q[0];
rz(-0.36488786) q[0];
rz(-0.94003135) q[1];
sx q[1];
rz(-0.53661984) q[1];
sx q[1];
rz(1.6392802) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4125) q[0];
sx q[0];
rz(-1.8209429) q[0];
sx q[0];
rz(-1.5604707) q[0];
x q[1];
rz(-0.83606007) q[2];
sx q[2];
rz(-1.9327455) q[2];
sx q[2];
rz(-1.8723633) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-2.9598436) q[1];
sx q[1];
rz(-2.1310924) q[1];
sx q[1];
rz(-2.3676141) q[1];
rz(3.0495166) q[3];
sx q[3];
rz(-1.3658938) q[3];
sx q[3];
rz(2.3931707) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-0.49729785) q[2];
sx q[2];
rz(-0.50540322) q[2];
sx q[2];
rz(-2.0765182) q[2];
rz(0.30125695) q[3];
sx q[3];
rz(-1.4091636) q[3];
sx q[3];
rz(1.6206954) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.168468) q[0];
sx q[0];
rz(-0.081806101) q[0];
sx q[0];
rz(0.43564963) q[0];
rz(-1.7565953) q[1];
sx q[1];
rz(-2.6711617) q[1];
sx q[1];
rz(-0.41697821) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0458841) q[0];
sx q[0];
rz(-2.4864712) q[0];
sx q[0];
rz(-1.6326293) q[0];
rz(-pi) q[1];
rz(2.4776393) q[2];
sx q[2];
rz(-1.10154) q[2];
sx q[2];
rz(-0.19367733) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(0.13874395) q[1];
sx q[1];
rz(-1.7094304) q[1];
sx q[1];
rz(2.0744051) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.74420332) q[3];
sx q[3];
rz(-2.9508698) q[3];
sx q[3];
rz(-2.423559) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-0.10432648) q[2];
sx q[2];
rz(-1.4929079) q[2];
sx q[2];
rz(-0.22582516) q[2];
rz(-0.2078235) q[3];
sx q[3];
rz(-2.4184629) q[3];
sx q[3];
rz(0.58661714) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
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
rz(2.726783) q[0];
sx q[0];
rz(-2.8746958) q[0];
sx q[0];
rz(1.6171932) q[0];
rz(2.1879451) q[1];
sx q[1];
rz(-1.8667659) q[1];
sx q[1];
rz(-1.3226002) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.246829) q[0];
sx q[0];
rz(-1.6920977) q[0];
sx q[0];
rz(1.1907207) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.0820504) q[2];
sx q[2];
rz(-2.6419123) q[2];
sx q[2];
rz(-2.3479455) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-0.10510124) q[1];
sx q[1];
rz(-2.0836012) q[1];
sx q[1];
rz(2.3829616) q[1];
rz(-pi) q[2];
rz(0.15354746) q[3];
sx q[3];
rz(-2.2205177) q[3];
sx q[3];
rz(-0.49680199) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.7913197) q[2];
sx q[2];
rz(-1.046317) q[2];
sx q[2];
rz(2.1255169) q[2];
rz(1.2223876) q[3];
sx q[3];
rz(-0.17761579) q[3];
sx q[3];
rz(2.5861752) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(pi/2) q[3];
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
rz(2.9983457) q[0];
sx q[0];
rz(-0.9220985) q[0];
sx q[0];
rz(-2.0621598) q[0];
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
