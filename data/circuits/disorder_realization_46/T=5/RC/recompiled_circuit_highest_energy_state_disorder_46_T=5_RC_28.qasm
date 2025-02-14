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
rz(-1.3396076) q[0];
sx q[0];
rz(-0.34914246) q[0];
sx q[0];
rz(-1.027663) q[0];
rz(-0.034962058) q[1];
sx q[1];
rz(5.265994) q[1];
sx q[1];
rz(7.979402) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5303237) q[0];
sx q[0];
rz(-1.9113386) q[0];
sx q[0];
rz(0.61113071) q[0];
x q[1];
rz(1.5367754) q[2];
sx q[2];
rz(-1.3933225) q[2];
sx q[2];
rz(2.3498775) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.6563182) q[1];
sx q[1];
rz(-2.7492011) q[1];
sx q[1];
rz(-0.57971422) q[1];
x q[2];
rz(-1.9361922) q[3];
sx q[3];
rz(-2.2789848) q[3];
sx q[3];
rz(-2.303249) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.33862904) q[2];
sx q[2];
rz(-0.45946071) q[2];
sx q[2];
rz(2.6098693) q[2];
rz(1.3945329) q[3];
sx q[3];
rz(-1.2732384) q[3];
sx q[3];
rz(-1.9664221) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8748473) q[0];
sx q[0];
rz(-1.6371472) q[0];
sx q[0];
rz(2.3496085) q[0];
rz(0.92164552) q[1];
sx q[1];
rz(-1.4896723) q[1];
sx q[1];
rz(-1.1080866) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4031206) q[0];
sx q[0];
rz(-2.5307025) q[0];
sx q[0];
rz(1.3752543) q[0];
rz(-pi) q[1];
x q[1];
rz(1.5032866) q[2];
sx q[2];
rz(-0.6246399) q[2];
sx q[2];
rz(2.4209674) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.0818357) q[1];
sx q[1];
rz(-1.1918157) q[1];
sx q[1];
rz(0.80767085) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.2543801) q[3];
sx q[3];
rz(-1.624148) q[3];
sx q[3];
rz(2.2406468) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(2.7532928) q[2];
sx q[2];
rz(-1.029195) q[2];
sx q[2];
rz(1.252906) q[2];
rz(-2.5168354) q[3];
sx q[3];
rz(-1.2478991) q[3];
sx q[3];
rz(-2.6784082) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6817634) q[0];
sx q[0];
rz(-2.9637931) q[0];
sx q[0];
rz(-0.27339992) q[0];
rz(-0.071648486) q[1];
sx q[1];
rz(-2.007808) q[1];
sx q[1];
rz(1.6260446) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8799389) q[0];
sx q[0];
rz(-2.0802705) q[0];
sx q[0];
rz(1.1306056) q[0];
x q[1];
rz(-0.63457527) q[2];
sx q[2];
rz(-0.56152841) q[2];
sx q[2];
rz(-0.47088366) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-1.3699941) q[1];
sx q[1];
rz(-2.150476) q[1];
sx q[1];
rz(-1.3808151) q[1];
x q[2];
rz(2.8876696) q[3];
sx q[3];
rz(-0.99231595) q[3];
sx q[3];
rz(1.8691269) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.7800954) q[2];
sx q[2];
rz(-2.4269035) q[2];
sx q[2];
rz(1.3698461) q[2];
rz(-2.0731549) q[3];
sx q[3];
rz(-1.3911824) q[3];
sx q[3];
rz(-0.95503241) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
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
rz(-0.84900981) q[0];
sx q[0];
rz(-2.7136901) q[0];
sx q[0];
rz(-0.12246116) q[0];
rz(-0.017008688) q[1];
sx q[1];
rz(-1.5127134) q[1];
sx q[1];
rz(2.2564127) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6017319) q[0];
sx q[0];
rz(-2.0575876) q[0];
sx q[0];
rz(2.9877325) q[0];
rz(2.3058168) q[2];
sx q[2];
rz(-1.3124497) q[2];
sx q[2];
rz(1.6222184) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.95751132) q[1];
sx q[1];
rz(-0.62465099) q[1];
sx q[1];
rz(2.2923325) q[1];
rz(0.2284352) q[3];
sx q[3];
rz(-0.71492787) q[3];
sx q[3];
rz(0.19107669) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(1.7077606) q[2];
sx q[2];
rz(-1.3154575) q[2];
sx q[2];
rz(-0.07746499) q[2];
rz(2.7205983) q[3];
sx q[3];
rz(-1.9546031) q[3];
sx q[3];
rz(0.25119701) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0780599) q[0];
sx q[0];
rz(-3.1227626) q[0];
sx q[0];
rz(-0.68191648) q[0];
rz(-2.6497427) q[1];
sx q[1];
rz(-2.1147155) q[1];
sx q[1];
rz(1.1600201) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.23719653) q[0];
sx q[0];
rz(-0.86441308) q[0];
sx q[0];
rz(-2.893413) q[0];
x q[1];
rz(3.0536041) q[2];
sx q[2];
rz(-2.5694429) q[2];
sx q[2];
rz(-2.3038626) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-0.3403389) q[1];
sx q[1];
rz(-1.7427708) q[1];
sx q[1];
rz(-2.1992963) q[1];
x q[2];
rz(-0.29472189) q[3];
sx q[3];
rz(-1.975276) q[3];
sx q[3];
rz(-1.3724788) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-2.457095) q[2];
sx q[2];
rz(-2.2726077) q[2];
sx q[2];
rz(-1.0082461) q[2];
rz(1.6393939) q[3];
sx q[3];
rz(-1.1049756) q[3];
sx q[3];
rz(-2.8777299) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.098792583) q[0];
sx q[0];
rz(-1.5941987) q[0];
sx q[0];
rz(2.7368326) q[0];
rz(-2.3233991) q[1];
sx q[1];
rz(-1.8408366) q[1];
sx q[1];
rz(2.4249446) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.19159) q[0];
sx q[0];
rz(-1.5609589) q[0];
sx q[0];
rz(-2.8623926) q[0];
rz(-2.0990853) q[2];
sx q[2];
rz(-2.2620438) q[2];
sx q[2];
rz(-1.3529568) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-2.4617052) q[1];
sx q[1];
rz(-2.6033834) q[1];
sx q[1];
rz(-0.52556441) q[1];
rz(-1.9394373) q[3];
sx q[3];
rz(-2.0799985) q[3];
sx q[3];
rz(-2.500755) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(1.3396259) q[2];
sx q[2];
rz(-2.5012987) q[2];
sx q[2];
rz(-2.4862945) q[2];
rz(-0.35704923) q[3];
sx q[3];
rz(-1.3909631) q[3];
sx q[3];
rz(2.6220139) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
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
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.94012564) q[0];
sx q[0];
rz(-0.31620142) q[0];
sx q[0];
rz(-2.1215718) q[0];
rz(0.54234281) q[1];
sx q[1];
rz(-2.0261363) q[1];
sx q[1];
rz(-1.382359) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4931204) q[0];
sx q[0];
rz(-0.97875226) q[0];
sx q[0];
rz(-0.73541321) q[0];
rz(-0.21127659) q[2];
sx q[2];
rz(-1.7045492) q[2];
sx q[2];
rz(1.1115896) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(0.51782896) q[1];
sx q[1];
rz(-1.6490855) q[1];
sx q[1];
rz(-2.0324367) q[1];
rz(3.0579733) q[3];
sx q[3];
rz(-1.573054) q[3];
sx q[3];
rz(1.797054) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(2.072784) q[2];
sx q[2];
rz(-1.6135608) q[2];
sx q[2];
rz(-2.7911348) q[2];
rz(-2.6751878) q[3];
sx q[3];
rz(-2.2240708) q[3];
sx q[3];
rz(2.1980227) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.35895178) q[0];
sx q[0];
rz(-0.39059165) q[0];
sx q[0];
rz(-2.1851831) q[0];
rz(1.6920754) q[1];
sx q[1];
rz(-2.3608975) q[1];
sx q[1];
rz(0.67109674) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.22204493) q[0];
sx q[0];
rz(-1.4853108) q[0];
sx q[0];
rz(1.7128904) q[0];
rz(-pi) q[1];
rz(2.8196149) q[2];
sx q[2];
rz(-0.59315943) q[2];
sx q[2];
rz(-2.4308487) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(2.8099212) q[1];
sx q[1];
rz(-1.5070032) q[1];
sx q[1];
rz(-1.924519) q[1];
rz(-1.3255886) q[3];
sx q[3];
rz(-0.88174654) q[3];
sx q[3];
rz(2.6264555) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(0.87390071) q[2];
sx q[2];
rz(-2.716422) q[2];
sx q[2];
rz(-1.1262061) q[2];
rz(2.8113484) q[3];
sx q[3];
rz(-1.4010022) q[3];
sx q[3];
rz(0.11708524) q[3];
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
rz(pi/2) q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0432334) q[0];
sx q[0];
rz(-2.5625304) q[0];
sx q[0];
rz(0.057057127) q[0];
rz(-0.13380274) q[1];
sx q[1];
rz(-1.0452784) q[1];
sx q[1];
rz(0.71279508) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0781435) q[0];
sx q[0];
rz(-2.8196555) q[0];
sx q[0];
rz(-2.9418895) q[0];
x q[1];
rz(2.7370896) q[2];
sx q[2];
rz(-1.1931915) q[2];
sx q[2];
rz(-1.1640358) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-2.5310933) q[1];
sx q[1];
rz(-1.9695909) q[1];
sx q[1];
rz(-0.51271768) q[1];
rz(-pi) q[2];
rz(-3.0607575) q[3];
sx q[3];
rz(-0.87817115) q[3];
sx q[3];
rz(2.5580542) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(1.628525) q[2];
sx q[2];
rz(-1.2863938) q[2];
sx q[2];
rz(-0.081550278) q[2];
rz(1.9214123) q[3];
sx q[3];
rz(-2.8704075) q[3];
sx q[3];
rz(-1.0903821) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.14209014) q[0];
sx q[0];
rz(-1.5220078) q[0];
sx q[0];
rz(-1.581544) q[0];
rz(-1.1355431) q[1];
sx q[1];
rz(-1.138843) q[1];
sx q[1];
rz(-1.9948237) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.97871507) q[0];
sx q[0];
rz(-0.93746569) q[0];
sx q[0];
rz(0.94453728) q[0];
x q[1];
rz(-2.3778474) q[2];
sx q[2];
rz(-2.1830705) q[2];
sx q[2];
rz(-2.6291922) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(0.15849491) q[1];
sx q[1];
rz(-1.6545516) q[1];
sx q[1];
rz(-0.2106481) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.0057529) q[3];
sx q[3];
rz(-1.6445789) q[3];
sx q[3];
rz(1.0467873) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-0.59209383) q[2];
sx q[2];
rz(-1.8213976) q[2];
sx q[2];
rz(2.7756694) q[2];
rz(-2.7538815) q[3];
sx q[3];
rz(-1.1185027) q[3];
sx q[3];
rz(1.4661192) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6112919) q[0];
sx q[0];
rz(-0.74321754) q[0];
sx q[0];
rz(2.8331953) q[0];
rz(-0.74956924) q[1];
sx q[1];
rz(-1.1309962) q[1];
sx q[1];
rz(1.8684734) q[1];
rz(-2.0012435) q[2];
sx q[2];
rz(-0.7067718) q[2];
sx q[2];
rz(1.5511647) q[2];
rz(-0.6324296) q[3];
sx q[3];
rz(-1.2865744) q[3];
sx q[3];
rz(-2.1435973) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
