OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(0.91166624) q[0];
sx q[0];
rz(-0.32719964) q[0];
sx q[0];
rz(1.3258452) q[0];
rz(1.3658547) q[1];
sx q[1];
rz(-0.39443016) q[1];
sx q[1];
rz(-2.3529513) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.74371007) q[0];
sx q[0];
rz(-1.8787483) q[0];
sx q[0];
rz(0.90713199) q[0];
rz(-pi) q[1];
rz(-2.2872988) q[2];
sx q[2];
rz(-0.88695705) q[2];
sx q[2];
rz(-1.8707866) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-0.44581283) q[1];
sx q[1];
rz(-1.7632293) q[1];
sx q[1];
rz(2.8532527) q[1];
rz(-0.12874916) q[3];
sx q[3];
rz(-2.1948619) q[3];
sx q[3];
rz(-0.20341132) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-0.49393645) q[2];
sx q[2];
rz(-1.6948573) q[2];
sx q[2];
rz(0.20230618) q[2];
rz(0.10937396) q[3];
sx q[3];
rz(-0.60905639) q[3];
sx q[3];
rz(3.0561395) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
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
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0789455) q[0];
sx q[0];
rz(-0.92342347) q[0];
sx q[0];
rz(1.1751291) q[0];
rz(1.4641948) q[1];
sx q[1];
rz(-2.1390095) q[1];
sx q[1];
rz(-2.7820803) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4660366) q[0];
sx q[0];
rz(-2.9630442) q[0];
sx q[0];
rz(1.9738865) q[0];
rz(-pi) q[1];
rz(0.56698842) q[2];
sx q[2];
rz(-0.63646454) q[2];
sx q[2];
rz(-1.3186962) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-3.0766616) q[1];
sx q[1];
rz(-0.97244579) q[1];
sx q[1];
rz(1.1538572) q[1];
rz(2.5595846) q[3];
sx q[3];
rz(-2.0872384) q[3];
sx q[3];
rz(-2.2565115) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-2.8466865) q[2];
sx q[2];
rz(-1.0470231) q[2];
sx q[2];
rz(2.8666551) q[2];
rz(1.8168195) q[3];
sx q[3];
rz(-1.5271657) q[3];
sx q[3];
rz(1.2980609) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0777968) q[0];
sx q[0];
rz(-0.15811385) q[0];
sx q[0];
rz(-2.7445444) q[0];
rz(-2.5322757) q[1];
sx q[1];
rz(-1.3343697) q[1];
sx q[1];
rz(1.5416001) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.9843053) q[0];
sx q[0];
rz(-1.2799036) q[0];
sx q[0];
rz(-0.10910927) q[0];
rz(1.7607401) q[2];
sx q[2];
rz(-3.0751588) q[2];
sx q[2];
rz(-0.22904598) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.065464547) q[1];
sx q[1];
rz(-2.4148126) q[1];
sx q[1];
rz(-1.5962334) q[1];
rz(2.2619369) q[3];
sx q[3];
rz(-2.1113951) q[3];
sx q[3];
rz(-2.4215655) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.93727532) q[2];
sx q[2];
rz(-1.5949564) q[2];
sx q[2];
rz(2.9079962) q[2];
rz(-1.8448081) q[3];
sx q[3];
rz(-1.9498884) q[3];
sx q[3];
rz(2.4209723) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.63491708) q[0];
sx q[0];
rz(-1.1577865) q[0];
sx q[0];
rz(1.7764212) q[0];
rz(-1.2719951) q[1];
sx q[1];
rz(-1.9422266) q[1];
sx q[1];
rz(-1.8796399) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.1256169) q[0];
sx q[0];
rz(-1.3787833) q[0];
sx q[0];
rz(-0.58108347) q[0];
x q[1];
rz(1.2620401) q[2];
sx q[2];
rz(-1.5456887) q[2];
sx q[2];
rz(-1.0756191) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(3.1168667) q[1];
sx q[1];
rz(-2.9477009) q[1];
sx q[1];
rz(0.34195447) q[1];
rz(-0.19660321) q[3];
sx q[3];
rz(-0.8340618) q[3];
sx q[3];
rz(1.7693258) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-2.9103553) q[2];
sx q[2];
rz(-2.1114025) q[2];
sx q[2];
rz(-0.78872284) q[2];
rz(-1.8606868) q[3];
sx q[3];
rz(-1.7773881) q[3];
sx q[3];
rz(-1.0219215) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
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
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7541517) q[0];
sx q[0];
rz(-1.0197637) q[0];
sx q[0];
rz(-2.4805241) q[0];
rz(-2.0263653) q[1];
sx q[1];
rz(-1.8293569) q[1];
sx q[1];
rz(0.95975319) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9929912) q[0];
sx q[0];
rz(-1.0558943) q[0];
sx q[0];
rz(1.881625) q[0];
x q[1];
rz(1.3044673) q[2];
sx q[2];
rz(-0.56364554) q[2];
sx q[2];
rz(0.061390419) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-0.57454039) q[1];
sx q[1];
rz(-2.8216719) q[1];
sx q[1];
rz(1.1362651) q[1];
rz(-2.2154551) q[3];
sx q[3];
rz(-1.7297267) q[3];
sx q[3];
rz(1.8862806) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(0.58012086) q[2];
sx q[2];
rz(-2.6802907) q[2];
sx q[2];
rz(2.0816154) q[2];
rz(2.3914242) q[3];
sx q[3];
rz(-1.3786022) q[3];
sx q[3];
rz(1.2257956) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.2724514) q[0];
sx q[0];
rz(-1.5831818) q[0];
sx q[0];
rz(0.22431746) q[0];
rz(-1.6250826) q[1];
sx q[1];
rz(-0.61619157) q[1];
sx q[1];
rz(0.075627653) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.62759864) q[0];
sx q[0];
rz(-2.450431) q[0];
sx q[0];
rz(-2.2014754) q[0];
x q[1];
rz(-1.2669771) q[2];
sx q[2];
rz(-1.4785789) q[2];
sx q[2];
rz(-0.95668787) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-1.5995248) q[1];
sx q[1];
rz(-2.0318527) q[1];
sx q[1];
rz(2.6510524) q[1];
rz(-pi) q[2];
x q[2];
rz(0.84556112) q[3];
sx q[3];
rz(-1.9048094) q[3];
sx q[3];
rz(-3.0507244) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(1.6077967) q[2];
sx q[2];
rz(-2.3888612) q[2];
sx q[2];
rz(3.094574) q[2];
rz(-2.4646344) q[3];
sx q[3];
rz(-1.2983863) q[3];
sx q[3];
rz(3.1162139) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.15636477) q[0];
sx q[0];
rz(-1.4893091) q[0];
sx q[0];
rz(1.697668) q[0];
rz(2.2185183) q[1];
sx q[1];
rz(-1.0052899) q[1];
sx q[1];
rz(2.9583171) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3175439) q[0];
sx q[0];
rz(-2.5185761) q[0];
sx q[0];
rz(-2.7035294) q[0];
rz(0.44281339) q[2];
sx q[2];
rz(-1.9573136) q[2];
sx q[2];
rz(-2.7903976) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-1.0456775) q[1];
sx q[1];
rz(-0.60904494) q[1];
sx q[1];
rz(-2.1629929) q[1];
rz(-pi) q[2];
rz(-2.8541862) q[3];
sx q[3];
rz(-1.0106405) q[3];
sx q[3];
rz(-1.0792102) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(3.1075333) q[2];
sx q[2];
rz(-1.6851765) q[2];
sx q[2];
rz(3.122186) q[2];
rz(0.16128811) q[3];
sx q[3];
rz(-1.015377) q[3];
sx q[3];
rz(1.2279145) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
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
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1041383) q[0];
sx q[0];
rz(-0.4087953) q[0];
sx q[0];
rz(-3.0492875) q[0];
rz(1.9654988) q[1];
sx q[1];
rz(-2.8832925) q[1];
sx q[1];
rz(-1.7324804) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.21299674) q[0];
sx q[0];
rz(-0.11830506) q[0];
sx q[0];
rz(-1.0932834) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.94495602) q[2];
sx q[2];
rz(-2.4279478) q[2];
sx q[2];
rz(-1.2882441) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.2947139) q[1];
sx q[1];
rz(-2.6637609) q[1];
sx q[1];
rz(0.61378693) q[1];
rz(0.61586611) q[3];
sx q[3];
rz(-1.0857673) q[3];
sx q[3];
rz(-2.6814493) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-3.1000503) q[2];
sx q[2];
rz(-1.5771733) q[2];
sx q[2];
rz(-1.1735631) q[2];
rz(-2.7581577) q[3];
sx q[3];
rz(-1.8337367) q[3];
sx q[3];
rz(1.0907382) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1293056) q[0];
sx q[0];
rz(-1.6467935) q[0];
sx q[0];
rz(-2.9360085) q[0];
rz(-2.3221817) q[1];
sx q[1];
rz(-1.2148379) q[1];
sx q[1];
rz(-0.49088556) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.67710971) q[0];
sx q[0];
rz(-0.77953952) q[0];
sx q[0];
rz(2.5960514) q[0];
x q[1];
rz(-2.728168) q[2];
sx q[2];
rz(-2.4707762) q[2];
sx q[2];
rz(2.6476423) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-3.1236022) q[1];
sx q[1];
rz(-1.523087) q[1];
sx q[1];
rz(-1.0022267) q[1];
rz(-pi) q[2];
x q[2];
rz(1.0596778) q[3];
sx q[3];
rz(-1.6569431) q[3];
sx q[3];
rz(2.3448157) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.72863355) q[2];
sx q[2];
rz(-1.0372936) q[2];
sx q[2];
rz(1.2566077) q[2];
rz(0.43577731) q[3];
sx q[3];
rz(-1.1083009) q[3];
sx q[3];
rz(0.88424879) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4202145) q[0];
sx q[0];
rz(-0.90737897) q[0];
sx q[0];
rz(0.23289982) q[0];
rz(-3.0026992) q[1];
sx q[1];
rz(-1.1537617) q[1];
sx q[1];
rz(0.5084261) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9121567) q[0];
sx q[0];
rz(-0.60750738) q[0];
sx q[0];
rz(-1.1077393) q[0];
rz(-pi) q[1];
x q[1];
rz(1.0632238) q[2];
sx q[2];
rz(-1.3295577) q[2];
sx q[2];
rz(-0.2080179) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-0.57181035) q[1];
sx q[1];
rz(-1.5628734) q[1];
sx q[1];
rz(1.7047764) q[1];
rz(-pi) q[2];
rz(-2.6026229) q[3];
sx q[3];
rz(-2.1930755) q[3];
sx q[3];
rz(-1.8405746) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-0.4109219) q[2];
sx q[2];
rz(-1.6795009) q[2];
sx q[2];
rz(2.3250568) q[2];
rz(-0.38980347) q[3];
sx q[3];
rz(-2.1392418) q[3];
sx q[3];
rz(-1.7635112) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6941866) q[0];
sx q[0];
rz(-0.92300713) q[0];
sx q[0];
rz(-0.22966455) q[0];
rz(1.5028839) q[1];
sx q[1];
rz(-1.8062183) q[1];
sx q[1];
rz(0.91581215) q[1];
rz(0.40789976) q[2];
sx q[2];
rz(-2.7277011) q[2];
sx q[2];
rz(-0.87521876) q[2];
rz(1.3744204) q[3];
sx q[3];
rz(-2.4300601) q[3];
sx q[3];
rz(-2.0236494) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
