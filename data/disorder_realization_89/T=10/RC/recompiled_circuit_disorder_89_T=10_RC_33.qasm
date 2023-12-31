OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(2.3857631) q[0];
sx q[0];
rz(-1.7321777) q[0];
sx q[0];
rz(0.29456079) q[0];
rz(0.28490588) q[1];
sx q[1];
rz(-0.51061881) q[1];
sx q[1];
rz(-2.7198305) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9298676) q[0];
sx q[0];
rz(-1.5063138) q[0];
sx q[0];
rz(0.026260016) q[0];
rz(-pi) q[1];
x q[1];
rz(2.4550081) q[2];
sx q[2];
rz(-1.6010487) q[2];
sx q[2];
rz(0.1498915) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-2.9383558) q[1];
sx q[1];
rz(-0.48523808) q[1];
sx q[1];
rz(2.7435702) q[1];
rz(-pi) q[2];
x q[2];
rz(0.59395091) q[3];
sx q[3];
rz(-1.6868601) q[3];
sx q[3];
rz(-0.56967294) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-1.1224147) q[2];
sx q[2];
rz(-0.77029595) q[2];
sx q[2];
rz(-2.1315234) q[2];
rz(-1.646237) q[3];
sx q[3];
rz(-1.3388747) q[3];
sx q[3];
rz(3*pi/11) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0533957) q[0];
sx q[0];
rz(-1.2258376) q[0];
sx q[0];
rz(-0.85900599) q[0];
rz(2.7711218) q[1];
sx q[1];
rz(-1.5444688) q[1];
sx q[1];
rz(-1.6765615) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5252936) q[0];
sx q[0];
rz(-1.787961) q[0];
sx q[0];
rz(-2.8263682) q[0];
x q[1];
rz(0.53497603) q[2];
sx q[2];
rz(-1.8638532) q[2];
sx q[2];
rz(1.6261634) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-0.49711455) q[1];
sx q[1];
rz(-1.2116355) q[1];
sx q[1];
rz(0.70980806) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.8584077) q[3];
sx q[3];
rz(-2.0809485) q[3];
sx q[3];
rz(0.19954296) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(1.7039965) q[2];
sx q[2];
rz(-1.2164755) q[2];
sx q[2];
rz(-0.55830467) q[2];
rz(1.0393418) q[3];
sx q[3];
rz(-1.2149518) q[3];
sx q[3];
rz(-2.5487652) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(-pi) q[3];
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
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.52477437) q[0];
sx q[0];
rz(-2.1891948) q[0];
sx q[0];
rz(2.9779789) q[0];
rz(2.4257461) q[1];
sx q[1];
rz(-0.84638458) q[1];
sx q[1];
rz(1.7680426) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.64840245) q[0];
sx q[0];
rz(-1.545854) q[0];
sx q[0];
rz(2.6148952) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.2839196) q[2];
sx q[2];
rz(-1.8048986) q[2];
sx q[2];
rz(0.26299325) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-0.92257567) q[1];
sx q[1];
rz(-0.95239988) q[1];
sx q[1];
rz(-2.4465357) q[1];
rz(-1.340629) q[3];
sx q[3];
rz(-1.7880511) q[3];
sx q[3];
rz(-2.1118856) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(0.13087656) q[2];
sx q[2];
rz(-2.5961582) q[2];
sx q[2];
rz(-2.1219357) q[2];
rz(0.9807469) q[3];
sx q[3];
rz(-0.5766944) q[3];
sx q[3];
rz(-2.5286634) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2402128) q[0];
sx q[0];
rz(-0.69832435) q[0];
sx q[0];
rz(-2.3024978) q[0];
rz(1.5377195) q[1];
sx q[1];
rz(-1.6042177) q[1];
sx q[1];
rz(-2.531321) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2903039) q[0];
sx q[0];
rz(-1.5820869) q[0];
sx q[0];
rz(1.4615061) q[0];
rz(-2.3158014) q[2];
sx q[2];
rz(-0.50160393) q[2];
sx q[2];
rz(0.85130168) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.8205386) q[1];
sx q[1];
rz(-2.39403) q[1];
sx q[1];
rz(1.1827521) q[1];
rz(-pi) q[2];
x q[2];
rz(2.2545492) q[3];
sx q[3];
rz(-2.4409557) q[3];
sx q[3];
rz(2.5168583) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.6976167) q[2];
sx q[2];
rz(-2.6306751) q[2];
sx q[2];
rz(-2.9042517) q[2];
rz(-0.90732968) q[3];
sx q[3];
rz(-1.5932339) q[3];
sx q[3];
rz(-1.2835519) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[3];
rz(-pi/2) q[3];
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
rz(-1.5089371) q[0];
sx q[0];
rz(-0.11015686) q[0];
sx q[0];
rz(-1.6089815) q[0];
rz(-0.73348796) q[1];
sx q[1];
rz(-1.2669022) q[1];
sx q[1];
rz(1.8283432) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8719296) q[0];
sx q[0];
rz(-2.4047244) q[0];
sx q[0];
rz(1.0760197) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.0079185) q[2];
sx q[2];
rz(-3.0745227) q[2];
sx q[2];
rz(1.9474533) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-0.65391958) q[1];
sx q[1];
rz(-0.54033414) q[1];
sx q[1];
rz(2.3565156) q[1];
rz(-pi) q[2];
rz(-0.58532183) q[3];
sx q[3];
rz(-0.96671852) q[3];
sx q[3];
rz(-0.6795336) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.0017172) q[2];
sx q[2];
rz(-1.2728609) q[2];
sx q[2];
rz(-2.4772947) q[2];
rz(0.70513606) q[3];
sx q[3];
rz(-2.4987529) q[3];
sx q[3];
rz(-3.0378708) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0649081) q[0];
sx q[0];
rz(-3.0397471) q[0];
sx q[0];
rz(-0.054811906) q[0];
rz(-2.1272155) q[1];
sx q[1];
rz(-1.6631815) q[1];
sx q[1];
rz(-0.48318133) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9478215) q[0];
sx q[0];
rz(-1.7521439) q[0];
sx q[0];
rz(0.062747196) q[0];
rz(-3.1088164) q[2];
sx q[2];
rz(-1.1083229) q[2];
sx q[2];
rz(1.0810766) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(2.3997765) q[1];
sx q[1];
rz(-0.92405926) q[1];
sx q[1];
rz(-0.6439376) q[1];
x q[2];
rz(-1.2383078) q[3];
sx q[3];
rz(-1.4629435) q[3];
sx q[3];
rz(-2.7319542) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-0.15423916) q[2];
sx q[2];
rz(-2.8581212) q[2];
sx q[2];
rz(-3.0563291) q[2];
rz(1.2003468) q[3];
sx q[3];
rz(-1.6253358) q[3];
sx q[3];
rz(2.9912662) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.36713704) q[0];
sx q[0];
rz(-3.0052003) q[0];
sx q[0];
rz(-2.1869587) q[0];
rz(0.57149354) q[1];
sx q[1];
rz(-1.1936455) q[1];
sx q[1];
rz(-2.8894997) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.68483401) q[0];
sx q[0];
rz(-1.7157468) q[0];
sx q[0];
rz(-0.42379871) q[0];
rz(-pi) q[1];
rz(0.7438296) q[2];
sx q[2];
rz(-0.93223909) q[2];
sx q[2];
rz(3.0628052) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-2.8512307) q[1];
sx q[1];
rz(-2.2154659) q[1];
sx q[1];
rz(2.6840997) q[1];
x q[2];
rz(-0.37766405) q[3];
sx q[3];
rz(-1.1389684) q[3];
sx q[3];
rz(-2.3218732) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(0.38348848) q[2];
sx q[2];
rz(-0.98758101) q[2];
sx q[2];
rz(-2.2371116) q[2];
rz(-2.1379743) q[3];
sx q[3];
rz(-2.2763054) q[3];
sx q[3];
rz(2.482567) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
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
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7082108) q[0];
sx q[0];
rz(-0.0032341783) q[0];
sx q[0];
rz(-1.7277539) q[0];
rz(0.66043234) q[1];
sx q[1];
rz(-1.3884037) q[1];
sx q[1];
rz(-1.7699014) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8577514) q[0];
sx q[0];
rz(-0.62566602) q[0];
sx q[0];
rz(1.7218504) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.76099446) q[2];
sx q[2];
rz(-1.4174263) q[2];
sx q[2];
rz(-2.209687) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.15118571) q[1];
sx q[1];
rz(-1.2745665) q[1];
sx q[1];
rz(-0.59466655) q[1];
rz(-pi) q[2];
x q[2];
rz(2.780704) q[3];
sx q[3];
rz(-1.4286255) q[3];
sx q[3];
rz(-1.7796381) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-2.3147605) q[2];
sx q[2];
rz(-2.5173352) q[2];
sx q[2];
rz(-0.39789847) q[2];
rz(2.6509616) q[3];
sx q[3];
rz(-1.5964973) q[3];
sx q[3];
rz(-1.3354966) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9234377) q[0];
sx q[0];
rz(-1.9240009) q[0];
sx q[0];
rz(-2.9072705) q[0];
rz(1.6802457) q[1];
sx q[1];
rz(-0.82273465) q[1];
sx q[1];
rz(-2.871002) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.71352495) q[0];
sx q[0];
rz(-1.3534091) q[0];
sx q[0];
rz(-2.9146951) q[0];
rz(-1.2241237) q[2];
sx q[2];
rz(-1.3895831) q[2];
sx q[2];
rz(-0.62435645) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-1.8447664) q[1];
sx q[1];
rz(-1.7080194) q[1];
sx q[1];
rz(-1.574135) q[1];
rz(1.7646592) q[3];
sx q[3];
rz(-2.6541775) q[3];
sx q[3];
rz(2.2262115) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.33971912) q[2];
sx q[2];
rz(-0.084771307) q[2];
sx q[2];
rz(1.6516997) q[2];
rz(0.70288944) q[3];
sx q[3];
rz(-1.1542164) q[3];
sx q[3];
rz(1.1606914) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7372195) q[0];
sx q[0];
rz(-1.5727366) q[0];
sx q[0];
rz(-2.9872966) q[0];
rz(-2.1986296) q[1];
sx q[1];
rz(-1.1318726) q[1];
sx q[1];
rz(-0.74238366) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.48001227) q[0];
sx q[0];
rz(-2.1614657) q[0];
sx q[0];
rz(2.3175879) q[0];
x q[1];
rz(-0.97787751) q[2];
sx q[2];
rz(-1.2061384) q[2];
sx q[2];
rz(-0.2702156) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-2.4780477) q[1];
sx q[1];
rz(-1.5596584) q[1];
sx q[1];
rz(-2.7958109) q[1];
rz(-pi) q[2];
rz(2.9418482) q[3];
sx q[3];
rz(-2.9962857) q[3];
sx q[3];
rz(-2.6655281) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-1.8538889) q[2];
sx q[2];
rz(-2.0159857) q[2];
sx q[2];
rz(-1.5819736) q[2];
rz(-2.93086) q[3];
sx q[3];
rz(-0.78453523) q[3];
sx q[3];
rz(2.6081086) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
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
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.29006526) q[0];
sx q[0];
rz(-2.1061438) q[0];
sx q[0];
rz(0.60594546) q[0];
rz(1.7998981) q[1];
sx q[1];
rz(-1.9683899) q[1];
sx q[1];
rz(1.2731332) q[1];
rz(2.8200091) q[2];
sx q[2];
rz(-0.84761534) q[2];
sx q[2];
rz(-1.9615704) q[2];
rz(2.3711754) q[3];
sx q[3];
rz(-2.9478248) q[3];
sx q[3];
rz(-0.41268681) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
